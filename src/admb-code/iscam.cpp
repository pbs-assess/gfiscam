	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#undef COUT
	#define COUT(object) cout << #object "\n" << object <<endl;
	#undef TINY
	#define TINY 1.e-08
	#undef NA
	#define NA -99.0
	#include <admodel.h>
	#include <time.h>
	#include <string.h>
	//#include "lib/ADMB_XMLDoc.h"
	#include "lib/msy.h"
	#include "lib/msy.hpp"
	#include "lib/baranov.h"
  #include "lib/LogisticNormal.h"
  #include "lib/milka.h"
  #include "lib/multinomial.h"
	#include "Selex.h"
	//#if defined _WIN32 || defined _WIN64
	#include "lib/LogisticNormal.cpp"
	#include "lib/LogisticStudentT.cpp"
	#include "lib/msy.cpp"
	#include "lib/baranov.cpp"
	#include "lib/milka.cpp"
	#include "lib/multinomial.cpp"
	#include "lib/multivariate_t.cpp"
	#include "lib/fvar_m58.cpp"
  //#endif
	ivector getIndex(const dvector& a, const dvector& b)
	{
		int i,j,n;
		n = 0;
		j = 1;
		for( i = a.indexmin(); i <= a.indexmax(); i++ )
		{
			 if(b(i) != NA) n++;
		}
		ivector tmp(1,n);
		for( i = a.indexmin(); i <= a.indexmax(); i++ )
		{
			if(b(i) != NA )
			{
				tmp(j++) = a(i);
			}
		}
		
		return(tmp);
	}
	void readMseInputs()
	{
	  	cout<<"yep this worked"<<endl;
	}
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	bool mcmcEvalPhase = 0;
	
	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;
	adstring stripExtension(adstring fileName)
	{
		/*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
		*/
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}
	
	
	// class selex_vector
	// {
	// 	private:
	// 		int m_length;
			
	// 	public:
	// 	selex_vector()
	// 	{
	// 		m_length = 0;
			
	// 	}
	// 	~selex_vector()
	// 	{
			
	// 	}
		
	// };
	
	// #ifdef __GNUDOS__
	//   #include <gccmanip.h>
	// #endif
	// Variables to store results from DIC calculations.
	double dicNoPar = 0;
	double dicValue = 0;
	  
	  
	
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <iscam.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  DataFile.allocate("DataFile");
  ControlFile.allocate("ControlFile");
  ProjectFileControl.allocate("ProjectFileControl");
  ProcedureControlFile.allocate("ProcedureControlFile");
  ScenarioControlFile.allocate("ScenarioControlFile");
 BaseFileName = stripExtension(ControlFile);  ///< BaseName given by the control file
 ReportFileName = BaseFileName + adstring(".rep");
 cout<<BaseFileName<<endl;
 ad_comm::change_datafile_name(ProjectFileControl);
  n_tac.allocate("n_tac");
 COUT(ProjectFileControl);
 COUT(n_tac);
  tac.allocate(1,n_tac,"tac");
  n_pfcntrl.allocate("n_pfcntrl");
  pf_cntrl.allocate(1,n_pfcntrl,"pf_cntrl");
  eof_pf.allocate("eof_pf");
		if(eof_pf!=-999)
		{
			cout<<"Error reading projection file."<<endl;
			cout<<"Last integer read is "<<eof_pf<<endl;
			cout<<"The file should end with -999.\n Aborting!"<<endl;
			ad_exit(1);
		}
		SimFlag=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retrieves the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			SimFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
		}
		
		// Catarina implementing a new command for generating new data control and pfc file
		// for a new project.
		NewFiles = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-new",opt))>-1)
		{
			NewFiles = 1;
			NewFileName = ad_comm::argv[on+1];
		}
		// command line option for retrospective analysis. "-retro retro_yrs"
		retro_yrs=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retro_yrs = atoi(ad_comm::argv[on+1]);
			cout<<"|____________________________________________________|\n";
			cout<<"| Implementing Retrospective analysis                |\n";
			cout<<"|____________________________________________________|\n";
			cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
		}
		// Management strategy evaluation.
		mseFlag = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt))>-1)
		{
			mseFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
			cout<<"|_________________________________________________|\n";
			cout<<"|Implementing Management Strategy Evaluation      |\n";
			cout<<"|_________________________________________________|\n";
		}
		// Test MSY
		testMSY = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-msy",opt))>-1)
		{
			cout<<"Testing MSY calculations with Spreadsheet MSF.xlsx"<<endl;
			testMSY = 1;
			
		}
 ad_comm::change_datafile_name(DataFile);
  narea.allocate("narea");
  ngroup.allocate("ngroup");
  nsex.allocate("nsex");
  syr.allocate("syr");
  nyr.allocate("nyr");
  sage.allocate("sage");
  nage.allocate("nage");
  ngear.allocate("ngear");
  age.allocate(sage,nage);
 n_ags = narea * ngroup * nsex;
 n_ag  = narea * ngroup;
 n_gs  = ngroup * nsex;
  n_area.allocate(1,n_ags);
  n_group.allocate(1,n_ags);
  n_sex.allocate(1,n_ags);
  pntr_ag.allocate(1,narea,1,ngroup);
  pntr_gs.allocate(1,ngroup,1,nsex);
  pntr_ags.allocate(1,narea,1,ngroup,1,nsex);
		age.fill_seqadd(sage,1);
		int ig,ih,is;
		ig = 0;
		ih = 0;
		is = 0;
		for(f=1; f<=narea; f++)
		{
			for(g=1; g<=ngroup; g++)
			{
				ih ++;
				pntr_ag(f,g) = ih;
				for(h=1;h<=nsex;h++)
				{
					ig ++;
					n_area(ig)  = f;
					n_group(ig) = g;
					n_sex(ig)   = h;
					pntr_ags(f,g,h) = ig;
					if(f==1)
					{
						is ++;
						pntr_gs(g,h) = is;
					}
				}
			}
		}
		if(!mseFlag)
		{
		cout<<"| ----------------------- |"<<endl;
		cout<<"| MODEL DIMENSION         |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| narea  \t"<<narea<<endl;
		cout<<"| ngroup \t"<<ngroup<<endl;
		cout<<"| nsex   \t"<<nsex<<endl;
		cout<<"| syr    \t"<<syr<<endl;
		cout<<"| nyr    \t"<<nyr<<endl;
		cout<<"| sage   \t"<<sage<<endl;
		cout<<"| nage   \t"<<nage<<endl;
		cout<<"| ngear  \t"<<ngear<<endl;
		cout<<"| n_area \t"<<n_area<<endl;
		cout<<"| n_group\t"<<n_group<<endl;
		cout<<"| n_sex  \t"<<n_sex<<endl;
		cout<<"| pntr_ag\n"<<pntr_ag<<endl;
		cout<<"| pntr_gs\n"<<pntr_gs<<endl;
		cout<<"| pntr_ags\n"<<pntr_ags(1)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		}
		
		/* Check for dimension errors in projection control file. */
		if( pf_cntrl(1)<syr || pf_cntrl(3)<syr || pf_cntrl(5)<syr )
		{
			cout<<"WARNING: start year in projection file control is" 
			" less than initial model year. Setting to syr."<<endl;
			// exit(1);
			pf_cntrl(1) = syr;
			pf_cntrl(3) = syr;
			pf_cntrl(5) = syr;
		}
		if( pf_cntrl(2)>nyr || pf_cntrl(4)>nyr || pf_cntrl(6)>nyr )
		{
			cout<<"ERROR: last year in projection file control is" 
			" greater than last model year."<<endl;
			exit(1);
		}
  dAllocation.allocate(1,ngear,"dAllocation");
  fsh_flag.allocate(1,ngear);
		dAllocation = dAllocation/sum(dAllocation);
		for(int k=1;k<=ngear;k++)
		{
			if(dAllocation(k)>0)
				fsh_flag(k)=1;
			else
				fsh_flag(k)=0;
		}
		nfleet = sum(fsh_flag);
  nFleetIndex.allocate(1,nfleet);
		j = 1;
		for(int k=1; k<=ngear;k++)
		{
			if(fsh_flag(k)) nFleetIndex(j++) = k;
		}
  d_linf.allocate(1,n_ags,"d_linf");
  d_vonbk.allocate(1,n_ags,"d_vonbk");
  d_to.allocate(1,n_ags,"d_to");
  d_a.allocate(1,n_ags,"d_a");
  d_b.allocate(1,n_ags,"d_b");
  d_ah.allocate(1,n_ags,"d_ah");
  d_gh.allocate(1,n_ags,"d_gh");
  n_MAT.allocate("n_MAT");
		if(n_MAT)
		{
			t1 = sage;
			t2 = nage;
		}
		else
		{
			t1 = 0;
			t2 = 0;
		}
  d_maturityVector.allocate(t1,t2,"d_maturityVector");
  la.allocate(1,n_ags,sage,nage);
  wa.allocate(1,n_ags,sage,nage);
  ma.allocate(1,n_ags,sage,nage);
		if(!mseFlag)
		{
		cout<<setw(8)<<setprecision(4)<<endl;
	  cout<<"| ----------------------- |"<<endl;
		cout<<"| GROWTH PARAMETERS       |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| d_linf  \t"<<d_linf<<endl;
	  	cout<<"| d_vonbk \t"<<d_vonbk<<endl;
	  	cout<<"| d_to    \t"<<d_to<<endl;
	  	cout<<"| d_a     \t"<<d_a<<endl;
	  	cout<<"| d_b     \t"<<d_b<<endl;
	  	cout<<"| d_ah    \t"<<d_ah<<endl;
	  	cout<<"| d_gh    \t"<<d_gh<<endl;
	  	cout<<"| ----------------------- |\n"<<endl;
	  	}
	  	// length & weight-at-age based on input growth pars
	  	ma.initialize();
                 for(ig=1;ig<=n_ags;ig++)
	  	{
	  		la(ig) = d_linf(ig)*(1. - exp(-d_vonbk(ig)*(age-d_to(ig))));
	  		wa(ig) = d_a(ig) * pow(la(ig),d_b(ig));
	  		h = n_sex(ig);
	  		if(n_MAT==0)
	  		{
	  			ma(ig) = plogis(age,d_ah(ig),d_gh(ig));
	  		}
	  		else if( n_MAT>0 && h !=2 )
	  		{
	  			ma(ig) = d_maturityVector;
	  		}
	  	}
		 
  nCtNobs.allocate("nCtNobs");
  dCatchData.allocate(1,nCtNobs,1,7,"dCatchData");
  d3_Ct.allocate(1,n_ags,syr,nyr,1,ngear);
		ft_count = nCtNobs;
		if(!mseFlag)
		{
		cout<<"| ----------------------- |"<<endl;
		cout<<"| HEAD(dCatchData)        |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<dCatchData.sub(1,3)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| TAIL(dCatchData)        |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<dCatchData.sub(nCtNobs-3,nCtNobs)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		}
		d3_Ct.initialize();
		int k;
		for(int ii=1;ii<=nCtNobs;ii++)
		{
			i = dCatchData(ii)(1);
			k = dCatchData(ii)(2);
			f = dCatchData(ii)(3);
			g = dCatchData(ii)(4);
			h = dCatchData(ii)(5);
			if( h==0 )
			{
				for(h=1;h<=nsex;h++)
				{
					ig = pntr_ags(f,g,h);
					d3_Ct(ig)(i)(k) = 1./nsex*dCatchData(ii)(7);
				}
			}
			else
			{
				ig = pntr_ags(f,g,h);
				d3_Ct(ig)(i)(k) = dCatchData(ii)(7);
			} 
		}
  nItNobs.allocate("nItNobs");
  nSurveyIndex.allocate(1,nItNobs);
  n_it_nobs.allocate(1,nItNobs,"n_it_nobs");
  n_survey_type.allocate(1,nItNobs,"n_survey_type");
  d3_survey_data.allocate(1,nItNobs,1,n_it_nobs,1,8,"d3_survey_data");
  it_wt.allocate(1,nItNobs,1,n_it_nobs);
  it_grp.allocate(1,nItNobs,1,n_it_nobs);
		if(!mseFlag)
		{
		cout<<"| ----------------------- |"<<endl;
		cout<<"| TAIL(d3_survey_data)       |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<d3_survey_data(nItNobs).sub(n_it_nobs(nItNobs)-3,n_it_nobs(nItNobs))<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		}
		for(int k=1;k<=nItNobs;k++)
		{
			it_wt(k) = column(d3_survey_data(k),7) + 1.e-30;
			it_grp(k)= column(d3_survey_data(k),5);
			nSurveyIndex(k) = d3_survey_data(k)(1,3);
		}
		double tmp_mu = mean(it_wt);
		for(int k=1;k<=nItNobs;k++)
		{
			it_wt(k) = it_wt(k)/tmp_mu;
		}
  nAgears.allocate("nAgears");
  n_A_nobs.allocate(1,nAgears,"n_A_nobs");
  n_A_sage.allocate(1,nAgears,"n_A_sage");
  n_A_nage.allocate(1,nAgears,"n_A_nage");
  inp_nscaler.allocate(1,nAgears,"inp_nscaler");
  n_ageFlag.allocate(1,nAgears,"n_ageFlag");
  d3_A.allocate(1,nAgears,1,n_A_nobs,n_A_sage-5,n_A_nage,"d3_A");
  d3_A_obs.allocate(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
		if( n_A_nobs(nAgears) > 0) 
		{
			if(!mseFlag)
			{
				if(n_A_nobs(nAgears) > 3)
				{
					cout<<"| ----------------------- |"<<endl;
					cout<<"| TAIL(A)       |"<<endl;
					cout<<"| ----------------------- |"<<endl;
					cout<<setw(4)<<d3_A(nAgears).sub(n_A_nobs(nAgears)-2,n_A_nobs(nAgears))<<endl;
					cout<<"| ----------------------- |\n"<<endl;
				  }
			}
			for(k=1;k<=nAgears;k++)
			{
				dmatrix tmp = trans(trans(d3_A(k)).sub(n_A_sage(k),n_A_nage(k)));
				if(inp_nscaler(k) > 0)
				{
					for(i = 1; i <= n_A_nobs(k); i++ )
					{
						 tmp(i) = tmp(i)/sum(tmp(i)) * inp_nscaler(k);
					}
				}
				d3_A_obs(k) = tmp;
				//d3_A_obs(k) = trans(trans(d3_A(k)).sub(n_A_sage(k),n_A_nage(k)));
				
			}      
		}
		else if(!mseFlag)
		{
			cout<<"| ----------------------- |"<<endl;
			cout<<"| NO AGE OR LENGTH DATA   |"<<endl;
			cout<<"| ----------------------- |"<<endl;
		}
  nWtTab.allocate("nWtTab");
  nWtNobs.allocate(1,nWtTab,"nWtNobs");
  d3_inp_wt_avg.allocate(1,nWtTab,1,nWtNobs,sage-5,nage,"d3_inp_wt_avg");
  tmp_nWtNobs.allocate(1,nWtTab);
  projwt.allocate(1,nWtTab);
  n_bf_wt_row.allocate(1,nWtTab);
  nMeanWt.allocate("nMeanWt");
  nMeanWtNobs.allocate(1,nMeanWt,"nMeanWtNobs");
  d3_mean_wt_data.allocate(1,nMeanWt,1,nMeanWtNobs,1,7,"d3_mean_wt_data");
		/*
		  This will determine the new dimension of d3_inp_wt_avg in case the backward 
		  projection is needed required and rename nWtNobs to tmp_nWtNobs 
		*/
		projwt.initialize();
		n_bf_wt_row.initialize();
		tmp_nWtNobs.initialize();
		for(int k=1; k<=nWtTab; k++)
		{
			tmp_nWtNobs(k) = nWtNobs(k);
			projwt(k)=1;
			for(i=1; i<=nWtNobs(k); i++)
			{
				if(nWtNobs(k) > 0 && d3_inp_wt_avg(k)(i)(sage-5) < 0)
				{
					n_bf_wt_row(k)++ ;
				}
			}	
			
			if(n_bf_wt_row(k)>0)
			{
				for(int i=1; i<=n_bf_wt_row(k); i++)
				{
					int exp_nyr = fabs(d3_inp_wt_avg(k,i,sage-5))-syr;
					tmp_nWtNobs(k) += exp_nyr; 
				}
				projwt(k)=-n_bf_wt_row(k);
			}	
				
			else if (n_bf_wt_row(k) == 0)
			{
				tmp_nWtNobs(k) = nWtNobs(k);
				projwt(k)=1;
			}
		}
		
		sum_tmp_nWtNobs = sum(tmp_nWtNobs);	
  xinp_wt_avg.allocate(1,nWtTab,1,tmp_nWtNobs,sage-5,nage);
  xxinp_wt_avg.allocate(1,sum_tmp_nWtNobs,sage-5,nage);
		/*
		  This will redimension the d3_inp_wt_avg  according to tmp_nWtNobs and rename 
		  the 3d array to xinp_wt_avg. Then the 3darray is converted to a matrix 
		  xxinp_wt_avg
		*/
		xinp_wt_avg.initialize();
		xxinp_wt_avg.initialize();
  		for(int k=1; k<=nWtTab; k++)
		{
			ivector iroww(0,n_bf_wt_row(k));
			iroww.initialize();
			if(nWtNobs(k) > 0)
			{
				if(n_bf_wt_row(k) > 0)
				{
					for(i=1; i<=n_bf_wt_row(k); i++)
					{
						d3_inp_wt_avg(k,i,sage-5) = fabs(d3_inp_wt_avg(k,i,sage-5));
						iroww(i) = d3_inp_wt_avg(k,i,sage-5)-syr+iroww(i-1);
						for(int jj=iroww(i);jj>=iroww(i-1)+1;jj--)
	 					{
	 						xinp_wt_avg(k)(jj)(sage-5) = syr+jj-iroww(i-1)-1 ;
	 						xinp_wt_avg(k)(jj)(sage-4,nage) = d3_inp_wt_avg(k)(i)(sage-4,nage);
	 					}
					}
					
					for(int jj = iroww(n_bf_wt_row(k))+1; jj <= tmp_nWtNobs(k); jj++)
	 				{
	 					xinp_wt_avg(k)(jj)(sage-5,nage) = d3_inp_wt_avg(k)(jj-iroww(n_bf_wt_row(k)))(sage-5,nage);
	 				}	
				}
				else
	 			{
	 				for(int jj = 1; jj <= tmp_nWtNobs(k); jj++)
	 				{
	 					xinp_wt_avg(k)(jj)(sage-5,nage) = d3_inp_wt_avg(k)(jj)(sage-5,nage);
	 				}
	 			}
		
				int ttmp =	sum(tmp_nWtNobs(1,k-1));
				int ttmp2 =	sum(tmp_nWtNobs(1,k));
				for(int jj=ttmp+1; jj<=ttmp2; jj++) 
				{
					xxinp_wt_avg(jj)(sage-5,nage) = xinp_wt_avg(k)(jj-ttmp)(sage-5,nage);
				}
			}
		}
  dWt_bar.allocate(1,n_ags,sage,nage);
  d3_wt_avg.allocate(1,n_ags,syr,nyr+1,sage,nage);
  d3_wt_dev.allocate(1,n_ags,syr,nyr+1,sage,nage);
  d3_wt_mat.allocate(1,n_ags,syr,nyr+1,sage,nage);
  d3_len_age.allocate(1,n_ags,syr,nyr+1,sage,nage);
		d3_wt_avg.initialize();
		d3_wt_dev.initialize();
		d3_wt_mat.initialize();
		d3_len_age.initialize();
		for(ig=1;ig<=n_ags;ig++)
		{
			for(int i = syr; i <= nyr; i++)
			{
				d3_wt_avg(ig)(i) = wa(ig);
				//d3_wt_mat(ig)(i) = pow(elem_prod(ma(ig),wa(ig)),d_iscamCntrl(6));
				d3_wt_mat(ig)(i) = elem_prod(ma(ig),wa(ig));
				d3_len_age(ig)(i) = pow(wa(ig)/d_a(ig),1./d_b(ig));
				// Insert calculations for ALK here.
			}
		}
		// the overwrite d3_wt_avg & d3_wt_mat with existing empirical data
		// SM Sept 6, 2013. Added option of using NA values (-99.0) for
		// missing weight-at-age data, or truncated age-data.
		int iyr;
		
		//if(nWtNobs(ii) > 0)
		//{
		for(i=1;i<=sum_tmp_nWtNobs;i++)
		{
			iyr = xxinp_wt_avg(i,sage-5);
			f   = xxinp_wt_avg(i,sage-3);
			g   = xxinp_wt_avg(i,sage-2);
			h   = xxinp_wt_avg(i,sage-1);
		// | SM Changed Sept 9, to accomodate NA's (-99) in empirical data.
			if( h )
			{
				ig                   = pntr_ags(f,g,h);
				dvector tmp          = xxinp_wt_avg(i)(sage,nage);
				ivector idx          = getIndex(age,tmp);
				for( int ii = 1; ii <= size_count(idx); ii++ )
				{
					d3_wt_avg(ig)(iyr)(idx(ii)) = xxinp_wt_avg(i)(idx(ii));
					d3_len_age(ig)(iyr)(idx(ii))= pow(d3_wt_avg(ig)(iyr)(idx(ii))
				                                  /d_a(ig),1./d_b(ig));
				 }
			//d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
			d3_wt_mat(ig)(iyr)      = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
			}
			else if( !h ) 
			{
				for(int h=1;h<=nsex;h++)
				{
					ig                   = pntr_ags(f,g,h);
					dvector tmp          = xxinp_wt_avg(i)(sage,nage);
					ivector idx          = getIndex(age,tmp);
					// Problem, array indexed differ, must loop over idx;
					// d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
					for( int ii = 1; ii <= size_count(idx); ii++)
					{
						d3_wt_avg(ig)(iyr)(idx(ii)) = xxinp_wt_avg(i)(idx(ii));
						d3_len_age(ig)(iyr)(idx(ii)) = pow(d3_wt_avg(ig)(iyr)(idx(ii))
					                               /d_a(ig),1./d_b(ig));
					}
					d3_wt_mat(ig)(iyr)      = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
				}
			}
		}
		//}
		
		// average weight-at-age in projection years
		for(ig=1;ig<=n_ags;ig++)
		{
			dWt_bar(ig)        = colsum(d3_wt_avg(ig).sub(pf_cntrl(3),pf_cntrl(4)));
			dWt_bar(ig)       /= pf_cntrl(4)-pf_cntrl(3)+1;
			d3_wt_avg(ig)(nyr+1) = dWt_bar(ig);
			d3_wt_mat(ig)(nyr+1) = elem_prod(dWt_bar(ig),ma(ig));
			d3_len_age(ig)(nyr+1) = pow(dWt_bar(ig)/d_a(ig),1./d_b(ig));
		}
		
		
		// deviations in mean weight-at-age
		for(ig=1;ig<=n_ags;ig++)
		{
			dmatrix mtmp = trans( d3_wt_avg(ig) );
			for(j=sage;j<=nage;j++)
			{
				//COUT(sum(first_difference(mtmp(j)(syr,nyr))));
				if( sum( first_difference(mtmp(j)(syr,nyr))) )
				{
					mtmp(j) = ( mtmp(j)-mean(mtmp(j)(syr,nyr)) ) 
							/ sqrt(var(mtmp(j)(syr,nyr)));
				}
				else
				{
					mtmp(j) = 0;
				}
			}
			d3_wt_dev(ig) = trans(mtmp);
			
			if( min(d3_wt_avg(ig))<=0.000 && min(d3_wt_avg(ig))!=NA )
			{
				cout<<"|-----------------------------------------------|"<<endl;
				cout<<"| ERROR IN INPUT DATA FILE FOR MEAN WEIGHT DATA |"<<endl;
				cout<<"|-----------------------------------------------|"<<endl;
				cout<<"| - Cannot have an observed mean weight-at-age  |"<<endl;
				cout<<"|   less than or equal to 0.  Please fix.       |"<<endl;
				cout<<"| - You are permitted to use '-99.0' for missing|"<<endl;
				cout<<"|   values in your weight-at-age data.          |"<<endl;
				cout<<"| - Aborting program!                           |"<<endl;
				cout<<"|-----------------------------------------------|"<<endl;
				ad_exit(1);
			}
		}
  eof.allocate("eof");
	  if(eof==999){
		cout<<"\n| -- END OF DATA SECTION -- |\n";
	  	cout<<"|         eof = "<<eof<<"         |"<<endl;
		cout<<"|___________________________|"<<endl;
	  }else{
		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	  }
  fmsy.allocate(1,ngroup,1,nfleet);
  fall.allocate(1,ngroup,1,nfleet);
  msy.allocate(1,ngroup,1,nfleet);
  bmsy.allocate(1,ngroup);
  age_tau2.allocate(1,nAgears);
 ad_comm::change_datafile_name(ControlFile);
  npar.allocate("npar");
  theta_control.allocate(1,npar,1,7,"theta_control");
  theta_ival.allocate(1,npar);
  theta_lb.allocate(1,npar);
  theta_ub.allocate(1,npar);
  theta_phz.allocate(1,npar);
  theta_prior.allocate(1,npar);
  ipar_vector.allocate(1,npar);
		theta_ival  = column(theta_control,1);
		theta_lb    = column(theta_control,2);
		theta_ub    = column(theta_control,3);
		theta_phz   = ivector(column(theta_control,4));
		theta_prior = ivector(column(theta_control,5));
		ipar_vector(1,2) = ngroup;
		ipar_vector(6,7) = ngroup;
		ipar_vector(3)   = n_gs;
		ipar_vector(4,5) = n_ag;
  nCompIndex.allocate(1,nAgears,"nCompIndex");
  nCompLikelihood.allocate(1,nAgears,"nCompLikelihood");
  dMinP.allocate(1,nAgears,"dMinP");
  dEps.allocate(1,nAgears,"dEps");
  nPhz_age_tau2.allocate(1,nAgears,"nPhz_age_tau2");
  nPhz_phi1.allocate(1,nAgears,"nPhz_phi1");
  nPhz_phi2.allocate(1,nAgears,"nPhz_phi2");
  nPhz_df.allocate(1,nAgears,"nPhz_df");
  check.allocate("check");
		if(check != -12345) 
		{
			COUT(check);cout<<"Error reading composition controls\n"<<endl; exit(1);
		}
  selex_controls.allocate(1,10,1,ngear,"selex_controls");
  isel_npar.allocate(1,ngear);
  jsel_npar.allocate(1,ngear);
  isel_type.allocate(1,ngear);
  sel_phz.allocate(1,ngear);
  n_sel_blocks.allocate(1,ngear);
  ahat_agemin.allocate(1,ngear);
  ghat_agemax.allocate(1,ngear);
  age_nodes.allocate(1,ngear);
  yr_nodes.allocate(1,ngear);
  lambda_1.allocate(1,ngear);
  lambda_2.allocate(1,ngear);
  lambda_3.allocate(1,ngear);
		ahat_agemin      = selex_controls(2);
		ghat_agemax      = selex_controls(3);
		age_nodes = selex_controls(4);
		yr_nodes  = selex_controls(5);
		lambda_1  = selex_controls(7);
		lambda_2  = selex_controls(8);
		lambda_3  = selex_controls(9);
		isel_type    = ivector(selex_controls(1));
		sel_phz      = ivector(selex_controls(6));
		n_sel_blocks = ivector(selex_controls(10));
  sel_blocks.allocate(1,ngear,1,n_sel_blocks,"sel_blocks");
		// | COUNT THE NUMBER OF ESTIMATED SELECTIVITY PARAMETERS TO ESTIMATE
		// | isel_npar number of columns for each gear.
		// | jsel_npar number of rows for each gear.
		isel_npar.initialize();
		for(i=1;i<=ngear;i++)
		{	
			jsel_npar(i)=1;
			switch(isel_type(i))
			{
				case 1:
					// logistic selectivity
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i); 
					break;
					
				case 2:
					// age-specific coefficients
					isel_npar(i) = (nage-sage);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 3:
				 	// cubic spline 
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 4:	 
					// annual cubic splines
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = (nyr-syr-retro_yrs)+1;
					break;
					
				case 5:  
					// bicubic spline
					jsel_npar(i) = age_nodes(i);
					isel_npar(i) = yr_nodes(i);
					break;
				
				case 6:
					// fixed logistic (no parameters estimated)
					// ensure sel_phz is set to negative value.
					isel_npar(i) = 2;
					if(sel_phz(i)>0) sel_phz(i) = -1;
					break;
					
				case 7:
					// CHANGED: Now working, Vivian Haist fixed it.
					// logistic (3 parameters) with mean body 
					// weight deviations. 
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 8:
					// Alternative logistic selectivity with d3_wt_dev coefficients.
					isel_npar(i) = 3;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 11:
					// Logistic length-based selectivity.
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 12:
					// Length-based selectivity coeffs with cubic spline interpolation
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = n_sel_blocks(i);
					break;
				case 13:
					// age-specific coefficients for agemin to agemax
					isel_npar(i) = (ghat_agemax(i)-ahat_agemin(i)+1);
					jsel_npar(i) = n_sel_blocks(i);
					
				default: break;
			}
		}
  nits.allocate("nits");
  q_prior.allocate(1,nits,"q_prior");
  mu_log_q.allocate(1,nits,"mu_log_q");
  sd_log_q.allocate(1,nits,"sd_log_q");
  fitMeanWt.allocate("fitMeanWt");
  nMeanWtCV.allocate("nMeanWtCV");
  weight_sig.allocate(1,nMeanWtCV,"weight_sig");
  d_iscamCntrl.allocate(1,15,"d_iscamCntrl");
  eofc.allocate("eofc");
		verbose = d_iscamCntrl(1);
		if(verbose) COUT(d_iscamCntrl);
		
		for(int ig=1;ig<=n_ags;ig++)
		{
			for(int i = syr; i <= nyr; i++)
			{
				d3_wt_mat(ig)(i) = pow(d3_wt_mat(ig)(i),d_iscamCntrl(6));
			}
		}
		if(eofc==999){
			cout<<"\n| -- END OF CONTROL SECTION -- |\n";
		  	cout<<"|          eofc = "<<eofc<<"          |"<<endl;
			cout<<"|______________________________|"<<endl;
		}else{
			cout<<"\n ***** ERROR CONTROL FILE ***** \n"<<endl; exit(1);
		}
  ilvec.allocate(1,8);
 ilvec    = ngear;
 ilvec(1) = 1;			
 ilvec(2) = nItNobs;			
 ilvec(3) = nAgears;		
 ilvec(4) = ngroup;
 ilvec(5) = nMeanWt;
  n_naa.allocate(1,nAgears);
 nyr = nyr - retro_yrs;
		if(retro_yrs)
		{
			if(pf_cntrl(2)>nyr) pf_cntrl(2) = nyr;
			if(pf_cntrl(4)>nyr) pf_cntrl(4) = nyr;
			if(pf_cntrl(6)>nyr) pf_cntrl(6) = nyr;
		}
		for( i = 1; i <= nCtNobs; i++ )
		{
			if( dCatchData(i)(1) > nyr ) ft_count --;
		}
		// Retrospective counter for n_A_nobs
		n_naa.initialize();
		for( k = 1; k <= nAgears; k++ )
		{
			for( i = 1; i <= n_A_nobs(k); i++ )
			{
				iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
				if( iyr <= nyr ) n_naa(k)++;
			}
		}
  n_saa.allocate(1,nAgears);
 syr = syr + (int)d_iscamCntrl(14);
	  for(int k = 1; k <= ngear; k++ ){
      sel_blocks(k)(1) = syr;
    }
    if(pf_cntrl(1)<syr) pf_cntrl(1) = syr;
    if(pf_cntrl(3)<syr) pf_cntrl(3) = syr;
    if(pf_cntrl(5)<syr) pf_cntrl(5) = syr;
    for( i = 1; i <= nCtNobs; i++ ){
      if( dCatchData(i)(1) < syr ) ft_count --;
    }
    // Prospective counter for n_A_nobs
    n_saa.initialize();
    n_saa = 1;
    for( int k = 1; k <= nAgears; k++ ){
      for( int i = 1; i <= n_A_nobs(k); i++ ){
        iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
        if( iyr < syr ){
          n_saa(k)++;
        }
      }
    }
 COUT((n_saa));
 COUT((n_naa));
 if(verbose) cout<<"||-- END OF DATA_SECTION --||"<<endl;
}

void model_parameters::initializationfunction(void)
{
  theta.set_initial_value(theta_ival);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  theta.allocate(1,npar,1,ipar_vector,theta_lb,theta_ub,theta_phz,"theta");
  sel_par.allocate(1,ngear,1,jsel_npar,1,isel_npar,-25.,25.,sel_phz,"sel_par");
		if ( !global_parfile )
		{
			for(int k=1; k<=ngear; k++)
			{
				if( isel_type(k)==1 || 
					isel_type(k)==6 || 
					(
					isel_type(k)>=7 && 
					isel_type(k) <= 12 
					)
					)
				{
					for(int j = 1; j <= n_sel_blocks(k); j++ )
					{
						double uu = 0;
						if(SimFlag && j > 1)
						{
							uu = 0.05*randn(j+rseed);
						} 
						sel_par(k,j,1) = log(ahat_agemin(k)*exp(uu));
						sel_par(k,j,2) = log(ghat_agemax(k));
					}
				}
				else if( isel_type(k) ==13 )
				{
					for(int j = 1; j <= n_sel_blocks(k); j++ )
					{
						double dd = 1.e-8;
						double stp = 1.0/(ghat_agemax(k)-ahat_agemin(k));
						sel_par(k)(j).fill_seqadd(dd,stp);
						//COUT(sel_par(k)(j));
						//exit(1);
					}
				}
			}
		}
  log_ft_pars.allocate(1,ft_count,-30.,3.0,1,"log_ft_pars");
		if(!SimFlag) log_ft_pars = log(0.10);
 int init_dev_phz = 2;
 if(d_iscamCntrl(5)) init_dev_phz = -1;
  init_log_rec_devs.allocate(1,n_ag,sage+1,nage,-15.,15.,init_dev_phz,"init_log_rec_devs");
  log_rec_devs.allocate(1,n_ag,syr,nyr,-15.,15.,2,"log_rec_devs");
 int m_dev_phz = -1;
     m_dev_phz = d_iscamCntrl(10);
 int  n_m_devs = d_iscamCntrl(12);
  log_m_nodes.allocate(1,n_m_devs,-5.0,5.0,m_dev_phz,"log_m_nodes");
  log_age_tau2.allocate(1,nAgears,-4.65,5.30,nPhz_age_tau2,"log_age_tau2");
  phi1.allocate(1,nAgears,-1.0,1.0,nPhz_phi1,"phi1");
  phi2.allocate(1,nAgears,0.0,1.0,nPhz_phi2,"phi2");
  log_degrees_of_freedom.allocate(1,nAgears,0.70,10.0,nPhz_df,"log_degrees_of_freedom");
  gamma_r.allocate(0,1,-4,"gamma_r");
gamma_r = 0;
  objfun.allocate("objfun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  m_bar.allocate("m_bar");
  #ifndef NO_AD_INITIALIZE
  m_bar.initialize();
  #endif
  ro.allocate(1,ngroup,"ro");
  #ifndef NO_AD_INITIALIZE
    ro.initialize();
  #endif
  bo.allocate(1,ngroup,"bo");
  #ifndef NO_AD_INITIALIZE
    bo.initialize();
  #endif
  sbo.allocate(1,ngroup,"sbo");
  #ifndef NO_AD_INITIALIZE
    sbo.initialize();
  #endif
  kappa.allocate(1,ngroup,"kappa");
  #ifndef NO_AD_INITIALIZE
    kappa.initialize();
  #endif
  steepness.allocate(1,ngroup,"steepness");
  #ifndef NO_AD_INITIALIZE
    steepness.initialize();
  #endif
  so.allocate(1,ngroup,"so");
  #ifndef NO_AD_INITIALIZE
    so.initialize();
  #endif
  beta.allocate(1,ngroup,"beta");
  #ifndef NO_AD_INITIALIZE
    beta.initialize();
  #endif
  m.allocate(1,n_gs,"m");
  #ifndef NO_AD_INITIALIZE
    m.initialize();
  #endif
  log_avgrec.allocate(1,n_ag,"log_avgrec");
  #ifndef NO_AD_INITIALIZE
    log_avgrec.initialize();
  #endif
  log_recinit.allocate(1,n_ag,"log_recinit");
  #ifndef NO_AD_INITIALIZE
    log_recinit.initialize();
  #endif
  q.allocate(1,nItNobs,"q");
  #ifndef NO_AD_INITIALIZE
    q.initialize();
  #endif
  ct.allocate(1,nCtNobs,"ct");
  #ifndef NO_AD_INITIALIZE
    ct.initialize();
  #endif
  eta.allocate(1,nCtNobs,"eta");
  #ifndef NO_AD_INITIALIZE
    eta.initialize();
  #endif
  log_m_devs.allocate(syr+1,nyr,"log_m_devs");
  #ifndef NO_AD_INITIALIZE
    log_m_devs.initialize();
  #endif
  rho.allocate(1,ngroup,"rho");
  #ifndef NO_AD_INITIALIZE
    rho.initialize();
  #endif
  varphi.allocate(1,ngroup,"varphi");
  #ifndef NO_AD_INITIALIZE
    varphi.initialize();
  #endif
  sig.allocate(1,ngroup,"sig");
  #ifndef NO_AD_INITIALIZE
    sig.initialize();
  #endif
  tau.allocate(1,ngroup,"tau");
  #ifndef NO_AD_INITIALIZE
    tau.initialize();
  #endif
  log_rt.allocate(1,n_ag,syr-nage+sage,nyr,"log_rt");
  #ifndef NO_AD_INITIALIZE
    log_rt.initialize();
  #endif
  nlvec.allocate(1,8,1,ilvec,"nlvec");
  #ifndef NO_AD_INITIALIZE
    nlvec.initialize();
  #endif
  epsilon.allocate(1,nItNobs,1,n_it_nobs,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  it_hat.allocate(1,nItNobs,1,n_it_nobs,"it_hat");
  #ifndef NO_AD_INITIALIZE
    it_hat.initialize();
  #endif
  qt.allocate(1,nItNobs,1,n_it_nobs,"qt");
  #ifndef NO_AD_INITIALIZE
    qt.initialize();
  #endif
  sbt.allocate(1,ngroup,syr,nyr+1,"sbt");
  #ifndef NO_AD_INITIALIZE
    sbt.initialize();
  #endif
  bt.allocate(1,ngroup,syr,nyr+1,"bt");
  #ifndef NO_AD_INITIALIZE
    bt.initialize();
  #endif
  rt.allocate(1,ngroup,syr+sage,nyr,"rt");
  #ifndef NO_AD_INITIALIZE
    rt.initialize();
  #endif
  delta.allocate(1,ngroup,syr+sage,nyr,"delta");
  #ifndef NO_AD_INITIALIZE
    delta.initialize();
  #endif
  annual_mean_weight.allocate(1,nMeanWt,1,nMeanWtNobs,"annual_mean_weight");
  #ifndef NO_AD_INITIALIZE
    annual_mean_weight.initialize();
  #endif
  obs_annual_mean_weight.allocate(1,nMeanWt,1,nMeanWtNobs,"obs_annual_mean_weight");
  #ifndef NO_AD_INITIALIZE
    obs_annual_mean_weight.initialize();
  #endif
  ft.allocate(1,n_ags,1,ngear,syr,nyr,"ft");
  #ifndef NO_AD_INITIALIZE
    ft.initialize();
  #endif
  F.allocate(1,n_ags,syr,nyr,sage,nage,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  M.allocate(1,n_ags,syr,nyr,sage,nage,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  Z.allocate(1,n_ags,syr,nyr,sage,nage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,n_ags,syr,nyr,sage,nage,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,n_ags,syr,nyr+1,sage,nage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  A_hat.allocate(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage,"A_hat");
  #ifndef NO_AD_INITIALIZE
    A_hat.initialize();
  #endif
  A_nu.allocate(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage,"A_nu");
  #ifndef NO_AD_INITIALIZE
    A_nu.initialize();
  #endif
  log_sel.allocate(1,ngear,1,n_ags,syr,nyr,sage,nage,"log_sel");
  #ifndef NO_AD_INITIALIZE
    log_sel.initialize();
  #endif
  sd_depletion.allocate(1,ngroup,"sd_depletion");
  sd_log_sbt.allocate(1,ngroup,syr,nyr+1,"sd_log_sbt");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
	// |---------------------------------------------------------------------------------|
	// | Run the model with input parameters to simulate real data.
	// |---------------------------------------------------------------------------------|
	// | - nf is a function evaluation counter.
 	// | - SimFlag comes from the -sim command line argument to simulate fake data.
 	// |
    nf=0;
  	if( testMSY )
  	{
  		testMSYxls();
  	}
	if( SimFlag ) 
	{
		initParameters();
		
		simulationModel(rseed);
	}
	
	if (NewFiles)
	{
		generate_new_files();	
	}
	
	if(verbose) cout<<"||-- END OF PRELIMINARY_CALCS_SECTION --||"<<endl;
	
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{100,  200,   500, 25000, 25000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{0.01, 0.01, 1.e-3, 1.e-4, 1.e-5}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
  objfun =0.0;
	initParameters();
	calcSelectivities(isel_type);
 	calcTotalMortality();
 	calcNumbersAtAge();
	calcTotalCatch();
	calcComposition();
	calcSurveyObservations();
	calcStockRecruitment();
	calcAnnualMeanWeight(); //START_RF_ADD ;  END_RF_ADD
	calcObjectiveFunction();
	if(sd_phase())
	{
		calcSdreportVariables();
	}
	if(mc_phase())
	{
		mcmcPhase=1;
	}
	if(mceval_phase())
	{
		mcmcEvalPhase=1;
		mcmc_output();
    cout<<"Running mceval phase"<<endl;
	}
	if( verbose ) {cout<<"End of main function calls"<<endl;}
}

void model_parameters::saveXMLFile(void)
{
	//ADMB_XMLDoc xml;
	/**
	Purpose:  This function calculates the sdreport variables.
	Author: Steven Martell
	Arguments:
		None
	NOTES:
	TODO list:
	  [?] - Calculate spawning biomass depletion for each group.
	*/
}

void model_parameters::calcSdreportVariables()
{
  {
	sd_depletion.initialize();
	sd_log_sbt.initialize();
	for(g=1;g<=ngroup;g++)
	{
		sd_depletion(g) = sbt(g)(nyr)/sbo(g);
		sd_log_sbt(g) = log(sbt(g));
	}
	if( verbose ) { cout<<"**** Ok after calcSdreportVariables ****"<<endl;}
  }
  	/**
  	Purpose: This function extracts the specific parameter values from the theta vector
  	       to initialize the leading parameters in the model.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  		- You must call this routine before running the simulation model to generate 
  		  fake data, otherwise you'll have goofy initial values for your leading parameters.
  		- Variance partitioning:
  	  Estimating total variance as = 1/precision
  	  and partition variance by rho = sig^2/(sig^2+tau^2).
  	  E.g. if sig = 0.2 and tau =1.12 then
  	  rho = 0.2^2/(0.2^2+1.12^2) = 0.03090235
  	  the total variance is kappa^2 = sig^2 + tau^2 = 1.2944
  	TODO list:
  	[ ] - Alternative parameterization using MSY and FMSY as leading parameters (Martell).
  	[*] - avg recruitment limited to area, may consider ragged object for area & stock.
  	*/
}

void model_parameters::initParameters()
{
  {
  	int ih;
	ro        = mfexp(theta(1));
	steepness = theta(2);
	m         = mfexp(theta(3));
	rho       = theta(6);
	varphi    = sqrt(1.0/theta(7));
	sig       = elem_prod(sqrt(rho) , varphi);
	tau       = elem_prod(sqrt(1.0-rho) , varphi);
	for(ih=1;ih<=n_ag;ih++)
	{
		log_avgrec(ih)  = theta(4,ih);
		log_recinit(ih) = theta(5,ih);
	}
	switch(int(d_iscamCntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			kappa = elem_div(4.*steepness,(1.-steepness));
			break;
		case 2:
			//Ricker model
			kappa = pow((5.*steepness),1.25);
		break;
	}
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
  }
}

dvar_vector model_parameters::cubic_spline(const dvar_vector& spline_coffs)
{
  {
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	//some testing here
	/*dvar_vector spline_nodes(1,nodes);
		spline_nodes.fill_seqadd(-0.5,1./(nodes-1));
		cout<<spline_nodes<<endl;
		vcubic_spline_function test_ffa(ia,spline_nodes);
		cout<<test_ffa(fa)<<endl;
		exit(1);*/
	return(ffa(fa));
  }
}

dvar_vector model_parameters::cubic_spline(const dvar_vector& spline_coffs, const dvector& la)
{
  {
	/*interplolation for length-based selectivity coefficeients*/
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(ffa(fa));
  }
  /**
   * @brief cubic spline interpolation
   * @details Uses cubic spline interpolatoin for data type variables based on a 
   * vector of spline coefficients, or nodes, and independent points.  
   * The nodes are rescaled to 0-1.  This function does not extrapolate beyond the 
   * independent points.
   * 
   * @param spline_coffs a data vector of spline coefficients (nodes)
   * @param la a vector of independent points for use in interpolation.
   * 
   * @return A data vector containing the interpolated points.
   */
}

dvector model_parameters::cubic_spline(const dvector& spline_coffs, const dvector& la)
{
  {
	/*interplolation for length-based selectivity coefficeients*/
	//RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	//RETURN_ARRAYS_DECREMENT();
	return(value(ffa(fa)));
	//return(1.0*la);
  }
  	/**
  	Purpose: This function loops over each of ngears and calculates the corresponding
  	         selectivity coefficients for that gear in each year.  It uses a switch 
  	         statement based on isel_type to determine which selectivty function to use
  	         for each particular gear that is specified in the control file.  See NOTES
  	         below for more information on selectivity models.
  	Author: Steven Martell
  	Arguments:
  		isel_type -> an ivector with integers that determine what selectivity model to use.
  	NOTES:
  		- The following is a list of the current selectivity models that are implemented:
		1)  Logistic selectivity with 2 parameters.
		2)  Age-specific selectivity coefficients with (nage-sage) parameters.
		    and the last two age-classes are assumed to have the same selectivity.
		3)  A reduced age-specific parameter set based on a bicubic spline.
		4)  Time varying cubic spline.
		5)  Time varying bicubic spline (2d version).
		6)  Fixed logistic.
		7)  Logistic selectivity based on relative changes in mean weight at age
		8)  Time varying selectivity based on logistic with deviations in 
		    weights at age (3 estimated parameters).
		11) Logistic selectivity with 2 parameters based on mean length.
		12) Length-based selectivity using cubic spline interpolation.
  		- The bicubic_spline function is located in stats.cxx library.
  	TODO list:
  	[*] add an option for length-based selectivity.  Use inverse of
		allometric relationship w_a = a*l_a^b; to get mean length-at-age from
		empirical weight-at-age data, then calculate selectivity based on 
		mean length. IMPLEMENTED IN CASE 11
	[*] change index for gear loop from j to k, and be consistent with year (i) and
	    age (j), and sex (h) indexing.
  	*/
}

void model_parameters::calcSelectivities(const ivector& isel_type)
{
  {
	int ig,i,j,k,byr,bpar,kgear;
	double tiny=1.e-10;
	dvariable p1,p2,p3;
	dvar_vector age_dev=age;
	dvar_matrix t1;
	dvar_matrix   tmp(syr,nyr-1,sage,nage);
	dvar_matrix  tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	// Selex cSelex(age);
	// logistic_selectivity cLogisticSelex(age);
	log_sel.initialize();
	for(kgear=1; kgear<=ngear; kgear++)
	{
		// The following is used to mirror another gear-type
		// based on the absolute value of sel_phz.
		k  = kgear;
		if(sel_phz(k) < 0)
		{
			k = abs(sel_phz(kgear));
			sel_par(kgear) = sel_par(k);
		}
		for( ig = 1; ig <= n_ags; ig++ )
		{
			tmp.initialize(); tmp2.initialize();
			dvector iy(1,yr_nodes(k));
			dvector ia(1,age_nodes(k));
			byr  = 1;
			bpar = 0; 
			switch(isel_type(k))
			{
				case 1: //logistic selectivity (2 parameters)
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						// cout<<"Testing selex class"<<endl;
						// log_sel(k)(ig)(i) = log( cSelex.logistic(sel_par(k)(bpar)) );
						// log_sel(k)(ig)(i) = log( cLogisticSelex(sel_par(k)(bpar)) );
						p1 = mfexp(sel_par(k,bpar,1));
						p2 = mfexp(sel_par(k,bpar,2));
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny );
					}
					break;
				case 6:	// fixed logistic selectivity
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					for(i=syr; i<=nyr; i++)
					{
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2) );
						// log_sel(k)(ig)(i) = log( cLogisticSelex(sel_par(k)(1)) );
					}
					break;
				case 2:	// age-specific selectivity coefficients
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						for(j=sage;j<=nage-1;j++)
						{
							log_sel(kgear)(ig)(i)(j)   = sel_par(k)(bpar)(j-sage+1);
						}
						log_sel(kgear)(ig)(i,nage) = log_sel(kgear)(ig)(i,nage-1);
					}
					break;
				case 3:	// cubic spline 
					for(i=syr; i<nyr; i++)
					{
						if( i==sel_blocks(k,byr) )
						{
							bpar ++;	
							log_sel(k)(ig)(i)=cubic_spline( sel_par(k)(bpar) );
							if( byr < n_sel_blocks(k) ) byr++;
						}
						log_sel(kgear)(ig)(i+1) = log_sel(k)(ig)(i);
					}
					break;
				case 4:	// time-varying cubic spline every year				
					for(i=syr; i<=nyr; i++)
					{
						log_sel(kgear)(ig)(i) = cubic_spline(sel_par(k)(i-syr+1));
					}
					break;
				case 5:	// time-varying bicubic spline
					ia.fill_seqadd( 0,1./(age_nodes(k)-1) );
					iy.fill_seqadd( 0,1./( yr_nodes(k)-1) );	
					bicubic_spline( iy,ia,sel_par(k),tmp2 );
					log_sel(kgear)(ig) = tmp2; 
					break;
				case 7:
					// time-varying selectivity based on deviations in weight-at-age
					// CHANGED This is not working and should not be used. (May 5, 2011)
					// SkDM:  I was not able to get this to run very well.
					// AUG 5, CHANGED so it no longer has the random walk component.
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					for(i = syr; i<=nyr; i++)
					{
						dvar_vector tmpwt=log(d3_wt_avg(ig)(i)*1000)/mean(log(d3_wt_avg(ig)*1000.));
						log_sel(kgear)(ig)(i) = log( plogis(tmpwt,p1,p2)+tiny );
					}	 
					break;
				case 8:
					//Alternative time-varying selectivity based on weight 
					//deviations (d3_wt_dev) d3_wt_dev is a matrix(syr,nyr+1,sage,nage)
					//p3 is the coefficient that describes variation in log_sel.
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					p3 = sel_par(k,1,3);
					for(i=syr; i<=nyr; i++)
					{
						tmp2(i) = p3*d3_wt_dev(ig)(i);
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny ) + tmp2(i);
					}
					break;
				case 11: // logistic selectivity based on mean length-at-age
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						p1 = mfexp(sel_par(k,bpar,1));
						p2 = mfexp(sel_par(k,bpar,2));
						dvector len = pow(d3_wt_avg(ig)(i)/d_a(ig),1./d_b(ig));
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(len,p1,p2) );
						//log_sel(kgear)(ig)(i) = log( plogis(len,p1,p2) );
					}	
					break;
				case 12: // cubic spline length-based coefficients.
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						dvector len = pow(d3_wt_avg(ig)(i)/d_a(ig),1./d_b(ig));
						log_sel(kgear)(ig)(i)=cubic_spline( sel_par(k)(bpar), len );
					}
					break;
					//parei aqui
				case 13:	// truncated age-specific selectivity coefficients
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						for(j=ahat_agemin(k); j<=ghat_agemax(k); j++)
						{
							log_sel(k)(ig)(i)(j)   = sel_par(k)(bpar)(j-ahat_agemin(k)+1);
						}
						for (j=ghat_agemax(k)+1; j<=nage; j++)
						{
							log_sel(kgear)(ig)(i,j) = log_sel(kgear)(ig)(i,ghat_agemax(k));
						}
						for(j=sage; j<ahat_agemin(k); j++)
						{
							log_sel(kgear)(ig)(i,j) = log_sel(kgear)(ig)(i,ahat_agemin(k));
						}						
					}
					break;
				default:
					log_sel(kgear)(ig)=0;
					break;
			}  // switch
			//subtract mean to ensure mean(exp(log_sel))==1
			//temp RF if statement
			//if(isel_type(k) !=6){
			//	for(i=syr;i<=nyr;i++)
			//	{
			//		log_sel(kgear)(ig)(i) -= log( mean(mfexp(log_sel(kgear)(ig)(i))) );
		       //
					// log_sel(k)(ig)(i) -= log( max(mfexp(log_sel(k)(ig)(i))) );
			//	}
			//} //end if	
		}
	}  //end of gear k
	if(verbose)cout<<"**** Ok after calcSelectivities ****"<<endl;
  }	
  	/**
  	Purpose: This function calculates fishing mortality, total mortality and annual
  	         surivival rates S=exp(-Z) for each age and year based on fishing mortality
  	         and selectivity coefficients.  Z also is updated with time-varying 
  	         natural mortality rates if specificed by user.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  		- Jan 5, 2012 Added catch_type to allow for catch in numbers, weight or spawn.
          In the case of spawn on kelp (roe fisheries), the Fishing mortality does not
          occur on the adult component.  
        - Added if(catch_type(k)!=3) //exclude roe fisheries
  		- F(group,year,age)
  		- Exclude type = 3, roe fisheries harvesting eggs only, not adults.
		- if dCatchData$sex is male & female combined, then allocate to both sexes.
  	TODO list:
  	[*] Dec 24, 2010.  Adding time-varying natural mortality.
  	[*] May 20, 2011.  Add cubic spline to the time-varying natural mortality.
	[ ] Calculate average M for reference point calculations based on pfc file.
	[ ] April 16, 2014. Adjust ft_count for retrospective and prospective analyses.
  	*/
}

void model_parameters::calcTotalMortality(void)
{
  {
	int ig,ii,i,k,l;
	int ft_counter = 0;
	dvariable ftmp;
	F.initialize(); 
	ft.initialize();
	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
       	for(ig=1;ig<=nCtNobs;ig++)
	{
		i  = dCatchData(ig)(1);	 //year
		k  = dCatchData(ig)(2);  //gear
		f  = dCatchData(ig)(3);  //area
		g  = dCatchData(ig)(4);  //group
		h  = dCatchData(ig)(5);  //sex
		l  = dCatchData(ig)(6);  //type
		if( i < syr ) continue;
		if( i > nyr ) continue;
		ft_counter ++;
		if( h )
		{
			ii = pntr_ags(f,g,h);    
			ftmp = mfexp(log_ft_pars(ft_counter));
			ft(ii)(k,i) = ftmp;
			if( l != 3 )
			{
				F(ii)(i) += ftmp*mfexp(log_sel(k)(ii)(i));
			}
		}
		else if( !h ) // h=0 case for asexual catch
		{
			for(h=1;h<=nsex;h++)
			{
				ii = pntr_ags(f,g,h);    
				ftmp = mfexp(log_ft_pars(ft_counter));
				ft(ii)(k,i) = ftmp;
				if( l != 3 )
				{
					F(ii)(i) += ftmp*mfexp(log_sel(k)(ii)(i));
				}		
			}
		}
	}
	// |---------------------------------------------------------------------------------|
	// | NATURAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// | - uses cubic spline to interpolate time-varying natural mortality
	M.initialize();
	log_m_devs.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		g = n_group(ig);
		h = n_sex(ig);
		M(ig) = m( pntr_gs(g,h) );
		if( active( log_m_nodes) )
		{
			int nodes = size_count(log_m_nodes);
			dvector im(1,nodes);
			dvector fm(syr+1,nyr);
			im.fill_seqadd(0,1./(nodes-1));
			fm.fill_seqadd(0,1./(nyr-syr));
			vcubic_spline_function m_spline(im,log_m_nodes);
			log_m_devs = m_spline( fm );
		}
		for(i=syr+1; i<=nyr; i++)
		{
			M(ig)(i) = M(ig)(i-1) * mfexp(log_m_devs(i));
		}
		// TODO fix for reference point calculations
		// m_bar = mean( M_tot.sub(pf_cntrl(1),pf_cntrl(2)) );	      
	}
	// |---------------------------------------------------------------------------------|
	// | TOTAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	  for(ig=1;ig<=n_ags;ig++)
	{
		Z(ig) = M(ig) + F(ig);
		S(ig) = mfexp(-Z(ig));
	}
	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;	
  }
  	/**
  	Purpose: This function initializes the numbers-at-age matrix in syr
  	         based on log_rinit and log_init_rec_devs, the annual recruitment
  	         based on log_rbar and log_rec_devs, and updates the number-at-age
  	         over time based on the survival rate calculated in calcTotalMortality.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
		- Aug 9, 2012.  Made a change here to initialize the numbers
		  at age in syr using the natural mortality rate at age in syr. 
		  Prior to this the average (m_bar) rate was used, since this 
		  has now changed with new projection control files.  Should only
		  affect models that were using time varying natural mortality.
  		- d_iscamCntrl(5) is a flag to start at unfished conditions, so set N(syr,sage) = ro
  	TODO list:
  	[ ] - Restrict log_avgrec and rec_devs to area and group dimensions (remove sex).
  	[ ] - Initialize from unfished conditions (d_iscamCntrl 5 flag is true then rt(syr) = ro)
  	*/
}

void model_parameters::calcNumbersAtAge(void)
{
  {
	int ig,ih;
	N.initialize();
	bt.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);
		dvar_vector lx(sage,nage);
		dvar_vector tr(sage,nage);
		lx(sage) = 1.0;
		for(j=sage;j< nage;j++)
		{
			lx(j+1) = lx(j) * exp( -M(ig)(syr)(j) );
		}
		lx(nage) /= (1.-exp(-M(ig)(syr,nage)));
		if( d_iscamCntrl(5) ) // initialize at unfished conditions.
		{
			tr =  log( ro(g) ) + log(lx);
		}
		else if ( !d_iscamCntrl(5) )
		{
			tr(sage)        = ( log_avgrec(ih)+log_rec_devs(ih)(syr));
			tr(sage+1,nage) = (log_recinit(ih)+init_log_rec_devs(ih));
			tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
		}
		N(ig)(syr)(sage,nage) = 1./nsex * mfexp(tr);
		log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);
		for(i=syr;i<=nyr;i++)
		{
			if( i>syr )
			{
				log_rt(ih)(i) = (log_avgrec(ih)+log_rec_devs(ih)(i));
				N(ig)(i,sage) = 1./nsex * mfexp( log_rt(ih)(i) );				
			}
			N(ig)(i+1)(sage+1,nage) =++elem_prod(N(ig)(i)(sage,nage-1)
			                                     ,S(ig)(i)(sage,nage-1));
			N(ig)(i+1,nage)        +=  N(ig)(i,nage)*S(ig)(i,nage);
			// average biomass for group in year i
			bt(g)(i) += N(ig)(i) * d3_wt_avg(ig)(i);
		}
		N(ig)(nyr+1,sage) = 1./nsex * mfexp( log_avgrec(ih));
		bt(g)(nyr+1) += N(ig)(nyr+1) * d3_wt_avg(ig)(nyr+1);
	}
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;	
  }
  	/**
  	Purpose:  This function calculates the predicted age-composition samples (A) for 
  	          both directed commercial fisheries and survey age-composition data. For 
  	          all years of data specified in the A matrix, calculated the predicted 
  	          proportions-at-age in the sampled catch-at-age.  If no catch-age data exist
  	          for a particular year i, for gear k (i.e. no directed fishery or from a 
  	          survey sample process that does not have an appreciable F), the calculate 
  	          the predicted proportion based on log(N) + log_sel(group,gear,year)
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  		- Adapted from iSCAM 1.5.  
  		- No longer using ragged arrays for gear, the ragged matrix is indexed by:
  		  year gear area, group, sex | age columns
  		- For the residuals, note that each gear is weigthed by the conditional MLE
  		  of the variance.
  	TODO list:
  	[x] - Merge redundant code from calcCatchAtAge
  	[*] - Add case where Chat data do not exsist.
	[x] - Calculate residuals A_nu; gets done automatically in dmvlogistic
	[?] - add plus group if n_A_nage < nage;  Aug 7, 2013
  	*/
}

void model_parameters::calcComposition(void)
{
  {
  	int ii,ig,kk;
  	dvar_vector va(sage,nage);
  	dvar_vector fa(sage,nage);
  	dvar_vector sa(sage,nage);
  	dvar_vector za(sage,nage);
  	dvar_vector ca(sage,nage);
  	dvar_vector na(sage,nage);
  	A_hat.initialize();
  	 for(kk=1;kk<=nAgears;kk++)
  	 {
  	 	for(ii=1;ii<=n_A_nobs(kk);ii++)
  	 	{
	  		i = d3_A(kk)(ii)(n_A_sage(kk)-5);
	  		k = d3_A(kk)(ii)(n_A_sage(kk)-4);
	  		f = d3_A(kk)(ii)(n_A_sage(kk)-3);
	  		g = d3_A(kk)(ii)(n_A_sage(kk)-2);
	  		h = d3_A(kk)(ii)(n_A_sage(kk)-1);
	  		// | trap for retrospecitve analysis.
	  		if(i < syr) continue;
	  		if(i > nyr) continue;
	  		if( h )  // age comps are sexed (h > 0)
	  		{
				ig = pntr_ags(f,g,h);
				va = mfexp(log_sel(k)(ig)(i));
				za = Z(ig)(i);
				sa = S(ig)(i);
				na = N(ig)(i);
				if( ft(ig)(k)(i)==0 )
				{
					ca = elem_prod(elem_prod(na,va),0.5*sa);
				}
				else
				{
					fa = ft(ig)(k)(i) * va;
					ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
				}
				//A_hat(kk)(ii) = ca(n_A_sage(kk),n_A_nage(kk));
				// | +group if n_A_nage(kk) < nage
				//if( n_A_nage(kk) < nage )
				//{
				//	A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
				//}
	  		}
	  		else if( !h )  // age-comps are unsexed
	  		{
	  			for(h=1;h<=nsex;h++)
	  			{
					ig = pntr_ags(f,g,h);
					va = mfexp(log_sel(k)(ig)(i));
					za = Z(ig)(i);
					sa = S(ig)(i);
					na = N(ig)(i);
					if( ft(ig)(k)(i)==0 )
					{
						ca = elem_prod(na,0.5*sa);
					}
					else
					{
						fa = ft(ig)(k)(i) * va;
						ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
					}
					//A_hat(kk)(ii) += ca(n_A_sage(kk),n_A_nage(kk));
					// | +group if n_A_nage(kk) < nage
					//if( n_A_nage(kk) < nage )
					//{
					//	A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
					//}
		  		}
	  		}
	  		// This is the age-composition
	  		if( n_ageFlag(kk) )
	  		{
	  			A_hat(kk)(ii) = ca(n_A_sage(kk),n_A_nage(kk));
	  			if( n_A_nage(kk) < nage )
				{
					A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
				}
	  		}
	  		else
	  		{
	  			/*
					This the catch-at-length composition.
					Pseudocode:
					-make an ALK
					-Ahat = ca * ALK
	  			*/
	  			dvar_vector mu = d3_len_age(ig)(i);
				dvar_vector sig= 0.1 * mu;
				dvector x(n_A_sage(kk),n_A_nage(kk));
				x.fill_seqadd(n_A_sage(kk),1);
				dvar_matrix alk = ALK(mu,sig,x);
	  			A_hat(kk)(ii) = ca * alk;
	  		}
	  		A_hat(kk)(ii) /= sum( A_hat(kk)(ii) );
  	 	}
  	}
	if(verbose)cout<<"**** Ok after calcComposition ****"<<endl;
  }	
}

void model_parameters::calcTotalCatch(void)
{
  {
  	/*
  	Purpose:  This function calculates the total catch.  
  	Dependencies: Must call calcCatchAtAge function first.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  	TODO list:
  	[ ] get rid of the obs_ct, ct, eta array structures, inefficient, better to use
  	    a matrix, then cbind the predicted catch and residuals for report. (ie. an R
  	    data.frame structure and use melt to ggplot for efficient plots.)
  	*/
  	int ii,l,ig;
  	double d_ct;
  	ct.initialize();
  	eta.initialize();
  	dvar_vector     fa(sage,nage);
  	dvar_vector     ca(sage,nage);
  	dvar_vector     sa(sage,nage);
  	dvar_vector     za(sage,nage);
  	for(ii=1;ii<=nCtNobs;ii++)
	{
		i    = dCatchData(ii,1);
		k    = dCatchData(ii,2);
		f    = dCatchData(ii,3);
		g    = dCatchData(ii,4);
		h    = dCatchData(ii,5);
		l    = dCatchData(ii,6);
		d_ct = dCatchData(ii,7);
  		// | trap for retro year
  		if( i<syr ) continue;
  		if( i>nyr ) continue;
		switch(l)
		{
			case 1:  // catch in weight
				if( h )
				{
					ig     = pntr_ags(f,g,h);
					fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
					za     = Z(ig)(i);
					sa     = S(ig)(i);
					ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
					ct(ii) = ca * d3_wt_avg(ig)(i);
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
						za     = Z(ig)(i);
						sa     = S(ig)(i);
						ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
						ct(ii)+= ca * d3_wt_avg(ig)(i);		
					}
				}
			break;
			case 2:  // catch in numbers
				if( h )
				{
					ig     = pntr_ags(f,g,h);
					fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
					za     = Z(ig)(i);
					sa     = S(ig)(i);
					ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
					ct(ii) = sum( ca );
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
						za     = Z(ig)(i);
						sa     = S(ig)(i);
						ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
						ct(ii)+= sum( ca );
					}
				}
			break;
			case 3:  // roe fisheries, special case
				if( h )
				{
					ig            = pntr_ags(f,g,h);
					dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
					ct(ii)        = (1.-exp(-ft(ig)(k)(i))) * ssb;
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig            = pntr_ags(f,g,h);
						dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
						ct(ii)       += (1.-exp(-ft(ig)(k)(i))) * ssb;
					}
				}
			break;
		}	// end of switch
		// | catch residual
		eta(ii) = log(d_ct+TINY) - log(ct(ii)+TINY);
	}
	if(verbose)cout<<"**** Ok after calcTotalCatch ****"<<endl;
  }
  	/**
  	Purpose:  This function computes the mle for survey q, calculates the survey 
  	          residuals (epsilon).
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  		- Oct 31, 2010, added retrospective counter.
  		- Nov 22, 2010, adding multiple surveys. 
  		Still need to check with retrospective option
  		- Nov 30, 2010, adjust the suvery biomass by the fraction of Z that has occurred 
  		when the survey was conducted. For herring spawning biomass this would be 
  		after the fishery has taken place.
  		- Dec 6, 2010, modified predicted survey biomass to accomodate empirical
  		weight-at-age data (d3_wt_avg).
  		- May 11, 2011.  Vivian Haist pointed out an error in survey biomass comparison.
  		The spawning biomass was not properly calculated in this routine. I.e. its 
  		different than the spawning biomass in the stock-recruitment routine. (Based on 
  		fecundity which changes with time when given empirical weight-at-age data.)
  		- Jan 6, 2012.  CHANGED corrected spawn survey observations to include a roe 
  		fishery that would remove potential spawn that would not be surveyed.
  		- d3_survey_data: (iyr index(it) gear area group sex wt timing)
  		- for MLE of survey q, using weighted mean of zt to calculate q.
  	TODO list:
  	    [?] - add capability to accomodate priors for survey q's.
  	    [ ] - verify q_prior=2 option for random walk in q.
  	    [ ] - For sel_type==3, may need to reduce abundance by F on spawning biomass (herring)
  	TODO LIST:
	  [ ] - add capability to accompodate priors for survey catchabiliyt coefficients.
  */
}

void model_parameters::calcSurveyObservations(void)
{
  {
	int ii,kk,ig,nz;
	double di;
	dvariable ftmp;
	dvar_vector Na(sage,nage);
	dvar_vector va(sage,nage);
	dvar_vector sa(sage,nage);
	epsilon.initialize();
	it_hat.initialize();
	for(kk=1;kk<=nItNobs;kk++)
	{
		// | Vulnerable number-at-age to survey.
		dvar_matrix V(1,n_it_nobs(kk),sage,nage);
		V.initialize();
		nz = 0;
		int iz=1;  // index for first year of data for prospective analysis.
		for(ii=1;ii<=n_it_nobs(kk);ii++)
		{
			i    = d3_survey_data(kk)(ii)(1);
			k    = d3_survey_data(kk)(ii)(3);
			f    = d3_survey_data(kk)(ii)(4);
			g    = d3_survey_data(kk)(ii)(5);
			h    = d3_survey_data(kk)(ii)(6);
			di   = d3_survey_data(kk)(ii)(8);
			// | trap for retrospective nyr change
			if( i < syr )
			{
				iz ++;
				nz ++;
				continue;
			} 
			if( i > nyr ) continue;
			nz ++;  // counter for number of observations.
			// h ==0?h=1:NULL;
			Na.initialize();
			for(h=1;h<=nsex;h++)
			{
				ig  = pntr_ags(f,g,h);
				va  = mfexp( log_sel(k)(ig)(i) );
				sa  = mfexp( -Z(ig)(i)*di );
				Na  = elem_prod(N(ig)(i),sa);
				switch(n_survey_type(kk))
				{
					case 1:
						V(ii) += elem_prod(Na,va);
					break; 
					case 2:
						V(ii) += elem_prod(elem_prod(Na,va),d3_wt_avg(ig)(i));
					break;
					case 3:
						V(ii) += elem_prod(Na,d3_wt_mat(ig)(i));
					break;
				}
			}
		} // end of ii loop
		dvector     it = trans(d3_survey_data(kk))(2)(iz,nz);
		dvector     wt = trans(d3_survey_data(kk))(7)(iz,nz);
		            wt = wt/sum(wt);
		dvar_vector t1 = rowsum(V);
		dvar_vector zt = log(it) - log(t1(iz,nz));
		dvariable zbar = sum(elem_prod(zt,wt));
				 q(kk) = mfexp(zbar);
		// | survey residuals
		epsilon(kk).sub(iz,nz) = zt - zbar;
		 it_hat(kk).sub(iz,nz) = q(kk) * t1(iz,nz);
		// | SPECIAL CASE: penalized random walk in q.
		if( q_prior(kk)==2 )
		{
			epsilon(kk).initialize();
			dvar_vector fd_zt     = first_difference(zt);
			dvariable  zw_bar     = sum(elem_prod(fd_zt,wt(iz,nz-1)));
			epsilon(kk).sub(iz,nz-1) = fd_zt - zw_bar;
			qt(kk)(iz) = exp(zt(iz));
			for(ii=iz+1;ii<=nz;ii++)
			{
				qt(kk)(ii) = qt(kk)(ii-1) * exp(fd_zt(ii-1));
			}
			it_hat(kk).sub(iz,nz) = elem_prod(qt(kk)(iz,nz),t1(iz,nz));
		}
	}
	if(verbose)cout<<"**** Ok after calcSurveyObservations ****"<<endl;
  }
  	/**
	Purpose:  
		This function is used to derive the underlying stock-recruitment 
		relationship that is ultimately used in determining MSY-based reference 
		points.  The objective of this function is to determine the appropriate 
		Ro, Bo and steepness values of either the Beverton-Holt or Ricker  Stock-
		Recruitment Model:
		Beverton-Holt Model
		\f$ Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau) \f$
		Ricker Model
		\f$ Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau) \f$
		The definition of a stock is based on group only. At this point, spawning biomass
		from all areas for a given group is the assumed stock, and the resulting
		recruitment is compared with the estimated recruits|group for all areas.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
		Psuedocode:
		-1) Get average natural mortality rate at age.
		-2) Calculate survivorship to time of spawning.
		-3) Calculate unfished spawning biomass per recruit.
		-4) Compute spawning biomass vector & substract roe fishery
		-5) Project spawning biomass to nyr+1 under natural mortality.
		-6) Calculate stock recruitment parameters (so, beta);
		-7) Calculate predicted recruitment
		-8) Compute residuals from estimated recruitments.
  	TODO list:
	  [?] - Change step 3 to be a weighted average of spawning biomass per recruit by area.
	  [?] - Increase dimensionality of ro, sbo, so, beta, and steepness to ngroup.
	  [?] - Add autocorrelation in recruitment residuals with parameter \f$ \gamma_r \f$.
  	*/
}

void model_parameters::calcStockRecruitment()
{
  {
  	int ig,ih;
  	rt.initialize();
  	sbt.initialize();
  	delta.initialize();
	dvariable phib;//,so,beta;
	dvector         fa(sage,nage); //fecundity here
	dvar_vector   stmp(sage,nage);
	dvar_vector     ma(sage,nage);
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector     lx(sage,nage); 
	dvar_vector     lw(sage,nage); 
	for(g=1;g<=ngroup;g++)
	{
		lx.initialize();
		lw.initialize();
		lx(sage) = 1.0;
		lw(sage) = 1.0;
		phib = 0;
		for(f=1;f<=narea;f++)
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);
				// | Step 1. average natural mortality rate at age.
				// | Step 2. calculate survivorship
				for(j=sage;j<=nage;j++)
				{
					ma(j) = mean(trans(M(ig))(j));
					fa(j) = mean( trans(d3_wt_mat(ig))(j) );
					if(j > sage)
					{
						lx(j) = lx(j-1) * mfexp(-ma(j-1));
					}
					lw(j) = lx(j) * mfexp(-ma(j)*d_iscamCntrl(13));
				}
				lx(nage) /= 1.0 - mfexp(-ma(nage));
				lw(nage) /= 1.0 - mfexp(-ma(nage));
				// | Step 3. calculate average spawing biomass per recruit.
				phib += 1./(narea*nsex) * lw*fa;
				// | Step 4. compute spawning biomass at time of spawning.
				for(i=syr;i<=nyr;i++)
				{
					stmp      = mfexp(-Z(ig)(i)*d_iscamCntrl(13));
					sbt(g,i) += elem_prod(N(ig)(i),d3_wt_mat(ig)(i)) * stmp;
				}
				// | Step 5. spawning biomass projection under natural mortality only.
				stmp          = mfexp(-M(ig)(nyr)*d_iscamCntrl(13));
				sbt(g,nyr+1) += elem_prod(N(ig)(nyr+1),d3_wt_mat(ig)(i)) * stmp;
			}
			// | Estimated recruits
			ih     = pntr_ag(f,g);
			rt(g) += mfexp(log_rt(ih)(syr+sage,nyr));
		}
		// | Step 6. calculate stock recruitment parameters (so, beta, sbo);
		so(g)  = kappa(g)/phib;
		sbo(g) = ro(g) * phib;
		// | Step 7. calculate predicted recruitment.
		dvar_vector tmp_st = sbt(g)(syr,nyr-sage).shift(syr+sage);
		switch(int(d_iscamCntrl(2)))
		{
			case 1:  // | Beverton Holt model
				beta(g)   = (kappa(g)-1.)/sbo(g);
				tmp_rt    = elem_div(so(g)*tmp_st,1.+beta(g)*tmp_st);
			break;
			case 2:  // | Ricker model
				beta(g)   = log(kappa(g))/sbo(g);
				tmp_rt    = elem_prod(so(g)*tmp_st,exp(-beta(g)*tmp_st));
			break;
		}
		// | Step 8. // residuals in stock-recruitment curve with gamma_r = 0
		delta(g) = log(rt(g))-log(tmp_rt)+0.5*tau(g)*tau(g);
		// Autocorrelation in recruitment residuals.
		// if gamma_r > 0 then 
		if( active(gamma_r) )
		{
			int byr = syr+sage+1;
			delta(g)(byr,nyr) 	= log(rt(g)(byr,nyr)) 
									- (1.0-gamma_r)*log(tmp_rt(byr,nyr)) 
									- gamma_r*log(++rt(g)(byr-1,nyr-1))
									+ 0.5*tau(g)*tau(g);			
		}
	}
     	if(verbose)cout<<"**** Ok after calcStockRecruitment ****"<<endl;
  }
  //RF added for comparison to Pacific Cod delay difference model	      //START_RF_ADD
}

void model_parameters::calcAnnualMeanWeight(void)
{
    {
    	     	int ii,kk,ig,nz;
		double di;
		dvar_vector wNa(sage,nage);
		dvar_vector wva(sage,nage);
		dvar_vector wsa(sage,nage);
		wNa.initialize();
		wva.initialize();
		wsa.initialize();
		for(kk=1;kk<=nMeanWt;kk++)   //loop through series with empirical annual mean weight data
		{
			dvar_matrix Vn(1,nMeanWtNobs(kk),sage,nage);	      // | Vulnerable number-at-age to gear
			dvar_matrix Vb(1,nMeanWtNobs(kk),sage,nage);	      // | Vulnerable biomass-at-age to gear
			Vn.initialize();
			Vb.initialize();
			nz = 0;
			int iz=1;  // index for first year of data for prospective analysis.
			for(ii=1;ii<=nMeanWtNobs(kk);ii++)	    //Loop through years 
			{
				i    = d3_mean_wt_data(kk)(ii)(1);  //year
				k    =d3_mean_wt_data(kk)(ii)(3);  //gear
				f    = d3_mean_wt_data(kk)(ii)(4);  //area
				g    = d3_mean_wt_data(kk)(ii)(5);  //group
				h    = d3_mean_wt_data(kk)(ii)(6);  //sex
				di   = d3_mean_wt_data(kk)(ii)(7); //timing
				// | trap for retrospective nyr change
				if( i < syr )
				{
					iz ++;
					nz ++;
					continue;
				} 
				if( i > nyr ) continue;
				nz ++;  // counter for number of observations.
				// h ==0?h=1:NULL;
				for(h=1;h<=nsex;h++)
				{
					ig  = pntr_ags(f,g,h);
					wva  = mfexp( log_sel(k)(ig)(i) );
					wsa  = mfexp( -Z(ig)(i)*di );   //accounts for survey timing
					wNa  = elem_prod(N(ig)(i),wsa);
					Vn(ii) += elem_prod(wNa,wva);  //adds sexes
					Vb(ii) += elem_prod(elem_prod(wNa,wva),d3_wt_avg(ig)(i));
				}
		annual_mean_weight(kk)(ii) = sum(Vb(ii))/sum(Vn(ii));
		obs_annual_mean_weight(kk)(ii)	= d3_mean_wt_data(kk)(ii)(2);	  //fill a matrix with observed annual mean weights - makes objective function calcs easier
		} // end of ii loop
	} // end of kk loop		
    	if(verbose)cout<<"**** Ok after calcAnnualMeanWeight ****"<<endl;
  } //end function	  //END_RF_ADD
}

void model_parameters::calcObjectiveFunction(void)
{
  {
  	/*
  	Purpose:  This function computes the objective function that ADMB will minimize.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
		There are several components to the objective function
		Likelihoods (nlvec):
			-1) likelihood of the catch data
			-2) likelihood of the survey abundance index
			-3) likelihood of age composition data 
			-4) likelihood for stock-recruitment relationship
			-5) penalized likelihood for fishery selectivities
			-6) penalized likelihood for fishery selectivities
			-7) penalized likelihood for fishery selectivities
			-8) likelihood for annual mean weight observations //START_RF_ADD   END_RF_ADD
  	TODO list:
	[*]	- Dec 20, 2010.  SJDM added prior to survey qs.
		  q_prior is an ivector with current options of 0 & 1 & 2.
		  0 is a uniform density (ignored) and 1 is a normal
		  prior density applied to log(q), and 2 is a random walk in q.
  	[ ] - Allow for annual sig_c values in catch data likelihood.
  	[ ] - Increase dimensionality of sig and tau to ngroup.
  	[ ] - Correct likelihood for cases when rho > 0 (Schnute & Richards, 1995)
  	*/
	nlvec.initialize();
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR CATCH DATA
	// |---------------------------------------------------------------------------------|
	// | - This likelihood changes between phases n-1 and n:
	// | - Phase (n-1): standard deviation in the catch based on user input d_iscamCntrl(3)
	// | - Phase (n)  : standard deviation in the catch based on user input d_iscamCntrl(4)
	// | 
	double sig_c =d_iscamCntrl(3);
	if(last_phase())
	{
		sig_c=d_iscamCntrl(4);
	}
	if( active(log_ft_pars) )
	{
		nlvec(1) = dnorm(eta,0.0,sig_c);
	}
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR RELATIVE ABUNDANCE INDICES
	// |---------------------------------------------------------------------------------|
	// | - sig_it     -> vector of standard deviations based on relative wt for survey.
	// |
	for(k=1;k<=nItNobs;k++)
	{
		ivector ig = it_grp(k);
		dvar_vector sig_it(1,n_it_nobs(k)); 
		for( i = 1; i <= n_it_nobs(k); i++ )
		{
			sig_it(i) = sig(ig(i))/it_wt(k,i);
		}
		nlvec(2,k)=dnorm(epsilon(k),sig_it);  
	}
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR AGE-COMPOSITION DATA
	// |---------------------------------------------------------------------------------|
	// | - Two options based on d_iscamCntrl(14):
	// | - 	1 -> multivariate logistic using conditional MLE of the variance for weight.
	// | -  2 -> multnomial, assumes input sample size as n in n log(p)
	// | -  Both likelihoods pool pmin (d_iscamCntrl(16)) into adjacent yearclass.
	// | -  PSEUDOCODE:
	// | -    => first determine appropriate dimensions for each of nAgears arrays (naa)
	// | -    => second extract sub arrays into obs (O) and predicted (P)
	// | -    => Compute either dmvlogistic, or dmultinom negative loglikehood.
	// | 
	// | TODO:
	// | [ ] - change A_nu to data-type variable, does not need to be differentiable.
	// | [ ] - issue 29. Fix submatrix O, P for prospective analysis & sex/area/group.
	A_nu.initialize();
	for(k=1;k<=nAgears;k++)
	{	
		if( n_A_nobs(k)>0 )
		{
			//int n_naa = 0;		//retrospective counter
			//int n_saa = 1;		//prospective counter
			int iyr;
			dmatrix      O(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
			dvar_matrix  P(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
			dvar_matrix nu(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
			O.initialize();
			P.initialize();
			nu.initialize();
			int ii=n_saa(k);
			for(i=1;i<=n_A_nobs(k);i++)
			{
				iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
				if(iyr >= syr && iyr <= nyr)
				{
					O(ii) = d3_A_obs(k)(i).sub(n_A_sage(k),n_A_nage(k));
					P(ii) = A_hat(k)(i).sub(n_A_sage(k),n_A_nage(k));
					ii ++;
				}
				//if( iyr <= nyr ) naa++;
				//if( iyr <  syr ) iaa++;
			}
			//dmatrix     O = trans(trans(d3_A_obs(k)).sub(n_A_sage(k),n_A_nage(k))).sub(iaa,naa);
			//dvar_matrix P = trans(trans(A_hat(k)).sub(n_A_sage(k),n_A_nage(k))).sub(iaa,naa);
			//dvar_matrix nu(O.rowmin(),O.rowmax(),O.colmin(),O.colmax()); 
			// | Choose form of the likelihood based on d_iscamCntrl(14) switch
			//switch(int(d_iscamCntrl(14)))
			logistic_normal cLN_Age( O,P,dMinP(k),dEps(k) );
			logistic_student_t cLST_Age( O,P,dMinP(k),dEps(k) );
			switch( int(nCompLikelihood(k)) )
			{
				case 1:
					nlvec(3,k) = dmvlogistic(O,P,nu,age_tau2(k),dMinP(k));
				break;
				case 2:
					nlvec(3,k) = dmultinom(O,P,nu,age_tau2(k),dMinP(k));
				break;
				case 3:
					if( !active(log_age_tau2(k)) )                 // LN1 Model
					{
						nlvec(3,k)  = cLN_Age();	
					}
					else
					{
						nlvec(3,k) = cLN_Age( exp(log_age_tau2(k)) );
					}
					// Residual
					if(last_phase())
					{
						nu          = cLN_Age.get_standardized_residuals();
						age_tau2(k) = cLN_Age.get_sigma2();
					}
				break;
				case 4:
					//logistic_normal cLN_Age( O,P,dMinP(k),dEps(k) );
					if( active(phi1(k)) && !active(phi2(k)) )  // LN2 Model
					{
            cout<<"\n\n\n\nLooking at log_age_tau2\n";
            cout<<log_age_tau2<<endl;
            cout<<"k="<<k<<", log_age_tau2(k)="<<log_age_tau2<<endl<<endl;
            cout<<"exp(log_age_tau2(k))="<<exp(log_age_tau2(k))<<endl<<endl;
            cout<<"phi2(k)="<<phi2(k)<<endl<<endl;
            cout<<"cLN_Age(expk,phi2k)="<<cLN_Age(exp(log_age_tau2(k)))<<endl<<endl;
						nlvec(3,k)   = cLN_Age(exp(log_age_tau2(k)),phi1(k));	
					}
					if( active(phi1(k)) && active(phi2(k)) )   // LN3 Model
					{
						nlvec(3,k)   = cLN_Age(exp(log_age_tau2(k)),phi1(k),phi2(k));	
					}
					// Residual
					if(last_phase())
					{
						nu          = cLN_Age.get_standardized_residuals();
						age_tau2(k) = cLN_Age.get_sigma2();
					}
				break;
				case 5: // Logistic-normal with student-t
					if( !active(log_degrees_of_freedom(k)) )
					{
						nlvec(3,k) = cLST_Age();
					}
					else
					{
						nlvec(3,k) = cLST_Age(exp(log_degrees_of_freedom(k)));
					}
					// Residual
					if(last_phase())
					{
						nu          = cLST_Age.get_standardized_residuals();
						age_tau2(k) = cLST_Age.get_sigma2();
					}
				break;
				case 6: // Multinomial with estimated effective sample size.
					nlvec(3,k) = mult_likelihood(O,P,nu,log_degrees_of_freedom(k));
				break; 
				case 7: // Multivariate-t 
					nlvec(3,k) = multivariate_t_likelihood(O,P,log_age_tau2(k),
					                                       log_degrees_of_freedom(k),
					                                       phi1(k),nu);
					age_tau2(k) = exp(value(log_age_tau2(k)));
				break;
			}
			// | Extract residuals.
			for(i=n_saa(k);i<=n_naa(k);i++)
			{
				A_nu(k)(i)(n_A_sage(k),n_A_nage(k))=nu(i);
			}
			ii = n_saa(k);
			for( i = 1; i <= n_A_nobs(k); i++ )
			{
				iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
				if(iyr >= syr && iyr <= nyr)
				{
					A_nu(k)(i)(n_A_sage(k),n_A_nage(k))=nu(ii++);		
				}
			}
		}
	}
	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT LIKELIHOOD COMPONENT
	// |---------------------------------------------------------------------------------|
	// | - tau is the process error standard deviation.
	if( active(theta(1)) || active(theta(2)) )
	{
		for(g=1;g<=ngroup;g++)
		{
			nlvec(4,g) = dnorm(delta(g),tau(g));
		}
	}
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD COMPONENTS FOR SELECTIVITY PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - lambda_1  -> penalty weight for smoothness
	// | - lambda_2  -> penalty weight for dome-shape
	// | - lambda_3  -> penalty weight for inter-annual variation.
	dvar_vector lvec(1,7); 
	lvec.initialize();
	int ig;
	for(k=1;k<=ngear;k++)
	{
		if(active(sel_par(k)))
		{
			//if not using logistic selectivity then
			//CHANGED from || to &&  May 18, 2011 Vivian
			if( isel_type(k)!=1 && 
				isel_type(k)!=7 && 
				isel_type(k)!=8 &&
				isel_type(k)!=11 )  
			{
				for(ig=1;ig<=n_ags;ig++)
				{
				for(i=syr;i<=nyr;i++)
				{
					//curvature in selectivity parameters
					dvar_vector df2 = first_difference(first_difference(log_sel(k)(ig)(i)));
					nlvec(5,k)     += lambda_1(k)/(nage-sage+1)*df2*df2;
					//penalty for dome-shapeness
					for(j=sage;j<=nage-1;j++)
						if(log_sel(k,ig,i,j)>log_sel(k,ig,i,j+1))
							nlvec(6,k)+=lambda_2(k)
										*square( log_sel(k,ig,i,j)-log_sel(k,ig,i,j+1) );
				}
				}
			}
			/*
			Oct 31, 2012 Halloween! Added 2nd difference penalty on time 
			for isel_type==(4)
			Mar 13, 2013, added 2nd difference penalty on isel_type==5 
			*/
			if( isel_type(k)==4 || isel_type(k)==5 || n_sel_blocks(k) > 1 )
			{
				for(ig=1;ig<=n_ags;ig++)
				{
				dvar_matrix trans_log_sel = trans( log_sel(k)(ig) );
				for(j=sage;j<=nage;j++)
				{
					dvar_vector df2 = first_difference(first_difference(trans_log_sel(j)));
					nlvec(7,k)     +=  lambda_3(k)/(nage-sage+1)*norm2(df2);
				}
				}
			}
		}
	}
	// |---------------------------------------------------------------------------------|
	// | CONSTRAINTS FOR SELECTIVITY DEVIATION VECTORS
	// |---------------------------------------------------------------------------------|
	// | [?] - TODO for isel_type==2 ensure mean 0 as well.
	// |
	for(k=1;k<=ngear;k++)
	{
		if( active(sel_par(k)) &&
			isel_type(k)!=1    &&
			isel_type(k)!=7    &&
			isel_type(k)!=8    &&
			isel_type(k)!=11 )
		{
			dvariable s = 0;
			if(isel_type(k)==5)  //bicubic spline version ensure column mean = 0
			{
				dvar_matrix tmp = trans(sel_par(k));
				for(j=1;j<=tmp.rowmax();j++)
				{
					s=mean(tmp(j));
					lvec(1)+=10000.0*s*s;
				}
			}
			if( isel_type(k)==2 ||
			    isel_type(k)==3 ||
			 	isel_type(k)==4 || 
				isel_type(k)==12 )
			{
				dvar_matrix tmp = sel_par(k);
				for(j=1;j<=tmp.rowmax();j++)
				{
					s=mean(tmp(j));
					lvec(1)+=10000.0*s*s;
				}
			}
		}
	}
	//START_RF_ADD
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR ANNUAL MEAN WEIGHT DATA
	// |---------------------------------------------------------------------------------|
	// | - sig_it     -> vector of standard deviations based on relative wt for survey.
	// |  init_3darray d3_mean_wt_data(1,nMeanWt,1,nMeanWtNobs,1,7);	
	for(k=1;k<=nMeanWt;k++)
	{
		dvar_vector epsilon_wt = log(annual_mean_weight(k)) - log(obs_annual_mean_weight(k));
		if(fitMeanWt) nlvec(8,k) = dnorm(epsilon_wt,weight_sig(k)); //fit to annual mean weight if fitMeanWt is switched on in the control file
	}
	//END_RF_ADD
	// |---------------------------------------------------------------------------------|
	// | PRIORS FOR LEADING PARAMETERS p(theta)
	// |---------------------------------------------------------------------------------|
	// | - theta_prior is a switch to determine which statistical distribution to use.
	// |
	dvariable ptmp; 
	dvar_vector priors(1,npar);
	priors.initialize();
	for(i=1;i<=npar;i++)
	{
		ptmp = 0;
		for(j=1;j<=ipar_vector(i);j++)
		{
			if( active(theta(i)) )
			{
				switch(theta_prior(i))
				{
				case 1:		//normal
					ptmp += dnorm(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
				case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
					ptmp += dlnorm(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
				case 3:		//beta distribution (0-1 scale)
					double lb,ub;
					lb=theta_lb(i);
					ub=theta_ub(i);
					ptmp += dbeta((theta(i,j)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
					break;
				case 4:		//gamma distribution
					ptmp += dgamma(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
				default:	//uniform density
					ptmp += log(1./(theta_control(i,3)-theta_control(i,2)));
					break;
				}
			}
		}
		priors(i) = ptmp;	
	}
	// |---------------------------------------------------------------------------------|
	// | PRIOR FOR SURVEY Q
	// |---------------------------------------------------------------------------------|
	// |
	dvar_vector qvec(1,nits);
	qvec.initialize();
	for(k=1;k<=nits;k++)
	{
		if(q_prior(k) == 1 )
		{
			qvec(k) = dnorm( log(q(k)), mu_log_q(k), sd_log_q(k) );
		}
	}
	// 	//** Legacy **  By accident took Rick Methot's bag from Nantes.
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD PENALTIES TO REGULARIZE SOLUTION
	// |---------------------------------------------------------------------------------|
	// | - pvec(1)  -> penalty on mean fishing mortality rate.
	// | - pvec(2)  -> penalty on first difference in natural mortality rate deviations.
	// | - pvec(4)  -> penalty on recruitment deviations.
	// | - pvec(5)  -> penalty on initial recruitment vector.
	// | - pvec(6)  -> constraint to ensure sum(log_rec_dev) = 0
	// | - pvec(7)  -> constraint to ensure sum(init_log_rec_dev) = 0
	// |
	dvar_vector pvec(1,7);
	pvec.initialize();
	dvariable log_fbar = mean(log_ft_pars);
	if(last_phase())
	{
		pvec(1) = dnorm(log_fbar,log(d_iscamCntrl(7)),d_iscamCntrl(9));
		// | Penalty for log_rec_devs (large variance here)
		for(g=1;g<=n_ag;g++)
		{
			pvec(4) += dnorm(log_rec_devs(g),2.0);
			pvec(5) += dnorm(init_log_rec_devs(g),2.0);
			dvariable s = 0;
			s = mean(log_rec_devs(g));
			pvec(6) += 1.e5 * s*s;
			s = mean(init_log_rec_devs(g));
			pvec(7) += 1.e5 * s*s;
		}
	}
	else
	{
		pvec(1) = dnorm(log_fbar,log(d_iscamCntrl(7)),d_iscamCntrl(8));
		//Penalty for log_rec_devs (CV ~ 0.0707) in early phases
		for(g=1;g<=n_ag;g++)
		{
			pvec(4) += 100.*norm2(log_rec_devs(g));
			pvec(5) += 100.*norm2(init_log_rec_devs(g));
			dvariable s = 0;
			s = mean(log_rec_devs(g));
			pvec(6) += 1.e5 * s*s;
			s = mean(init_log_rec_devs(g));
			pvec(7) += 1.e5 * s*s;
		}
	}
	if(active(log_m_nodes))
	{
		double std_mdev = d_iscamCntrl(11);
		dvar_vector fd_mdevs=first_difference(log_m_devs);
		pvec(2)  = dnorm(fd_mdevs,std_mdev);
		pvec(2) += 0.5*norm2(log_m_nodes);
	}
	if(verbose)
	{
		COUT(nlvec);
		COUT(lvec);
		COUT(priors);
		COUT(pvec);
		COUT(qvec);
	}
	// COUT(nlvec);
	objfun  = sum(nlvec);
	objfun += sum(lvec);
	objfun += sum(priors);
	objfun += sum(pvec);
	objfun += sum(qvec);
	nf++;
	if(verbose)cout<<"**** Ok after calcObjectiveFunction ****"<<endl;
	//cout<<"nlvec3 is "<<nlvec(3)<<endl;
	//	ad_exit(1); //para!
  }
}

void model_parameters::calcReferencePoints()
{
  {
  	/*
  	Purpose:  This function calculates the MSY-based reference points, and also loops
  	          over values of F and F-multipliers and calculates the equilibrium yield
  	          for each fishing gear.
  	Author: Steven Martell
  	Arguments:
  		None
  	NOTES:
  		- This function is based on the msyReferencePoint class object written by 
  		  Steve Martell on the Island of Maui while on Sabbatical leave from UBC.
  		- The msyReferencePoint class uses Newton-Raphson method to iteratively solve
  		  for the Fmsy values that maximize catch for each fleet. You can compare 
  		  MSY-reference points with the numerical values calculated at the end of this
  		  function.
  		- Code check: appears to find the correct value of MSY
		  in terms of maximizing ye.  Check to ensure rec-devs
		  need a bias correction term to get this right.
		- Modification for multiple fleets:
    	  	Need to pass a weighted average vector of selectivities
    	  	to the equilibrium routine, where the weights for each
    	  	selectivity is based on the dAllocation to each fleet.
	 		Perhaps as a default, assign an equal dAllocation to each
	 		fleet.  Eventually,user must specify dAllocation in 
	 		control file.  DONE
	 	- Use selectivity in the terminal year to calculate reference
	 	  points.  See todo, this is something that should be specified by the user.
		- June 8, 2012.  SJDM.  Made the following changes to this routine.
			1) changed reference points calculations to use the average
			   weight-at-age and fecundity-at-age.
			2) change equilibrium calculations to use the catch dAllocation
			   for multiple gear types. Not the average vulnerablity... this was wrong.
		- July 29, 2012.  SJDM Issue1.  New routine for calculating reference points
		  for multiple fleets. In this case, finds a vector of Fmsy's that simultaneously 
		  maximizes the total catch for each of the fleets respectively.  See
		  iSCAMequil_soln.R for an example.
		- August 1, 2012.  SJDM, In response to Issue1. A major overhaul of this routine.
		  Now using the new Msy class to calculate reference points. This greatly simplifies
		  the code in this routine and makes other routines (equilibrium) redundant.  Also
		  the new Msy class does a much better job in the case of multiple fleets.
		- Aug 11, 2012.
		  For the Pacific herring branch omit the get_fmsy calculations and use only the 
		  Bo calculuation for the reference points.  As there are no MSY based reference
		  points required for the descision table. 
		- May 8, 2013.
		  Starting to modify this code to allow for multiple areas, sex and stock-specific
		  reference points.  Key to this modification is the definition of a stock. For
		  the purposes of reference points, a stock is assumed to be distributed over all
		  areas, can be unisex or two sex, and can be fished by all fleets given the stock
		  exists in an area where the fishing gear operates.
		- Aug, 3-6, 2013.
		  Major effort went into revising this routine as well as the Msy class to 
		  calculate MSY-based reference points for multiple fleets and multiple sexes.  
		  I had found a significant bug in the dye calculation where I need to use the 
		  proper linear algebra to calculate the vector of dye values. This appears to 
		  be working properly and I've commented out the lines of code where I numerically
		  checked the derivatives of the catch equation.  This is a major acomplishment.
		- Mar, 2013.
		  A major new development here with the use of msy.hpp and a template class for 
		  calculating MSY-based reference points.  The user can now calculate reference
		  points for each gear based on fixed allocation, and optimum allocations based
		  on relative differences in selectivities among the gears landing fish. Uses the
		  name space "rfp".
   	PSEUDOCODE: 
   		(1) : Construct array of selectivities (potentially sex based log_sel)
   		(2) : Construct arrays of d3_wt_avg and d3_wt_mat for reference years.
	  	(3) : Come up with a reasonable guess for fmsy for each gear in nfleet.
	  	(4) : Instantiate an Msy class object and get_fmsy.
	  	(5) : Use Msy object to get reference points.
	slx::Selex<dvar_vector> * ptr;  //Pointer to Selex base class
  ptr = new slx::LogisticCurve<dvar_vector,dvariable>(mu,sd);
  log_sel = ptr->logSelectivity(age);
  delete ptr;
  	TODO list:
  	[ ] - allow user to specify which selectivity years are used in reference point
  	      calculations. This should probably be done in the projection File Control.
  	*/
	int kk,ig;
	// | (1) : Matrix of selectivities for directed fisheries.
	// |     : log_sel(gear)(n_ags)(year)(age)
	// |     : ensure dAllocation sums to 1.
	dvector d_ak(1,nfleet);
	d3_array  d_V(1,n_ags,1,nfleet,sage,nage);
	dvar3_array  dvar_V(1,n_ags,1,nfleet,sage,nage);
	for(k=1;k<=nfleet;k++)
	{
		kk      = nFleetIndex(k);
		d_ak(k) = dAllocation(kk);
		for(ig=1;ig<=n_ags;ig++)
		{
			d_V(ig)(k) = value( exp(log_sel(kk)(ig)(nyr)) );
			dvar_V(ig)(k) =( exp(log_sel(kk)(ig)(nyr)) );
		}
	}
	d_ak /= sum(d_ak);
	// | (2) : Average weight and mature spawning biomass for reference years
	// |     : dWt_bar(1,n_ags,sage,nage)
	dmatrix fa_bar(1,n_ags,sage,nage);
	dmatrix  M_bar(1,n_ags,sage,nage);
	for(ig=1;ig<=n_ags;ig++)
	{
		fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
		M_bar(ig)  = colsum(value(M(ig).sub(pf_cntrl(3),pf_cntrl(4))));
		M_bar(ig) /= pf_cntrl(4)-pf_cntrl(3)+1;	
	}
	// | (3) : Initial guess for fmsy for each fleet
	// |     : set fmsy = 2/3 of M divided by the number of fleets
	fmsy.initialize();
	fall.initialize();
	msy.initialize();
	bmsy.initialize();
	dvar_vector dftry(1,nfleet);
	dftry  = 0.6/nfleet * mean(M_bar);
	// | (4) : Instantiate msy class for each stock
	for(g=1;g<=ngroup;g++)
	{
		double d_rho = d_iscamCntrl(13);
		dvector d_mbar = M_bar(g);
		dvector   d_wa = dWt_bar(g);
		dvector   d_fa = fa_bar(g);
		//Pointer to the base class
		//rfp::referencePoints<dvariable,dvar_vector,dvar_matrix> * pMSY; 
		//pMSY = new rfp::msy<dvariable,dvar_vector,dvar_matrix,dvar3_array>
		//(ro(g),steepness(g),d_rho,M_bar,dWt_bar,fa_bar,dvar_V);
		//dvar_vector dfmsy = pMSY->getFmsy(dftry);
		//delete pMSY;
		//cout<<"Initial Fe "<<dftry<<endl;		 //RF turned this off
		rfp::msy<dvariable,dvar_vector,dvar_matrix,dvar3_array> 
		c_MSY(ro(g),steepness(g),d_rho,M_bar,dWt_bar,fa_bar,dvar_V);
		dvar_vector dfmsy = c_MSY.getFmsy(dftry,d_ak);
		bo  = c_MSY.getBo();
		dvariable dbmsy = c_MSY.getBmsy();
		dvar_vector dmsy = c_MSY.getMsy();
		bmsy(g) = value(dbmsy);
		msy(g)  = value(dmsy);
		fmsy(g) = value(dfmsy);
		//c_MSY.print();	    //RF turned this off
	}
	// Data-type version of MSY-based reference points.
	for( ig = 1; ig <= n_ags; ig++ )
	{
		fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
		M_bar(ig)  = colsum(value(M(ig).sub(pf_cntrl(3),pf_cntrl(4))));
		M_bar(ig) /= pf_cntrl(4)-pf_cntrl(3)+1;	
	}
	for( g = 1; g <= ngroup; g++ )
	{
		double d_ro  = value(ro(g));
		double d_h   = value(steepness(g));
		double d_rho = d_iscamCntrl(13);
		rfp::msy<double,dvector,dmatrix,d3_array>
		c_dMSY(d_ro,d_h,d_rho,M_bar,dWt_bar,fa_bar,d_V);
		fmsy(g) = c_dMSY.getFmsy(value(dftry));
		bo = c_dMSY.getBo();
		bmsy(g) = c_dMSY.getBmsy();
		msy(g)  = c_dMSY.getMsy();
		//c_dMSY.print();    RF turned this off
		dvector finit(1,nfleet);
		finit=fmsy(g);
		c_dMSY.checkDerivatives(finit);
		//cout<<"group \t"<<g<<endl;
		//exit(1);
		Msy c_msy(d_ro,d_h,M_bar,d_rho,dWt_bar,fa_bar,&d_V);
		fmsy(g) = 0.1;
		c_msy.get_fmsy(fmsy(g));
		bo = c_msy.getBo();
		bmsy(g) = c_msy.getBmsy();
		msy(g) = c_msy.getMsy();
		//cout<<"Old Msy class"<<endl;	   RF turned this off
		//c_msy.print();	  RF turned this off
	}
	if(verbose)cout<<"**** Ok after calcReferencePoints ****"<<endl;
  }
  /**
   * This is a simple test routine for comparing the MSY class output to the 
   * MSF.xlsx spreadsheet that was used to develop the multiple fleet msy 
   * method.  Its a permanent feature of the iscam code for testing.
   */	
}

void model_parameters::testMSYxls()
{
	double ro = 1.0;
	double steepness = 0.75;
	double d_rho = 0.0;
	dmatrix m_bar(1,1,1,20);
	m_bar = 0.30;
	dmatrix dWt_bar(1,1,1,20);
	dWt_bar(1).fill("{0.005956243,0.035832542,0.091848839,0.166984708,0.252580458,0.341247502,0.427643719,0.508367311,0.581557922,0.646462315,0.703062177,0.751788795,0.793319224,0.828438536,0.857951642,0.882630331,0.903184345,0.920248191,0.934377843,0.946053327}");
	dmatrix fa_bar(1,1,1,20);
	fa_bar(1).fill("{0,0,0,0,0.252580458,0.341247502,0.427643719,0.508367311,0.581557922,0.646462315,0.703062177,0.751788795,0.793319224,0.828438536,0.857951642,0.882630331,0.903184345,0.920248191,0.934377843,0.946053327}");
	d3_array dvar_V(1,1,1,2,1,20);
	dvar_V(1)(1).fill("{0.001271016,0.034445196,0.5,0.965554804,0.998728984,0.999954602,0.99999838,0.999999942,0.999999998,1,1,1,1,1,1,1,1,1,1,1}");
	dvar_V(1)(2).fill("{0.000552779,0.006692851,0.07585818,0.5,0.92414182,0.993307149,0.999447221,0.999954602,0.999996273,0.999999694,0.999999975,0.999999998,1,1,1,1,1,1,1,1}");
	//dvar_V(1)(2).fill("{0.000189406,0.000789866,0.003287661,0.013576917,0.054313266,0.19332137,0.5,0.80667863,0.945686734,0.986423083,0.996712339,0.999210134,0.999810594,0.999954602,0.99998912,0.999997393,0.999999375,0.99999985,0.999999964,0.999999991}");
	dvector dftry(1,2);
	dftry = 0.1 ;
	cout<<"Initial Fe "<<dftry<<endl;	    
	rfp::msy<double,dvector,dmatrix,d3_array> 
	c_MSY(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dfmsy = c_MSY.getFmsy(dftry);
	cout<<"Fmsy = "<<dfmsy<<endl;		  
	dvector ak(1,2);
	ak = 0.3;
	ak(2) = 1-ak(1);
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSYk(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dkmsy = c_MSYk.getFmsy(dftry,ak);
	cout<<"Fmsy_k ="<<dkmsy<<endl;
	c_MSYk.print();
	dvector akmsy = c_MSYk.getFmsy(dftry);
	c_MSYk.print();
	exit(1);
  	/**
  	Purpose:  This routine gets called from the PRELIMINARY_CALCS_SECTION if the 
  	          user has specified the -sim command line option.  The random seed
  	          is specifed at the command line.
  	Author: Steven Martell
  	Arguments:
  		seed -> a random seed for generating a unique, repeatable, sequence of 
  		        random numbers to be used as observation and process errors.
  	NOTES:
		- This routine will over-write the observations in memory
		  with simulated data, where the true parameter values are
		  the initial values.  Change the standard deviations of the 
		  random number vectors epsilon (observation error) or 
 		  recruitment devs wt (process error).
 		- At the end of DATA_SECTION nyrs is modified by retro_yrs if -retro.
 		- Add back the retro_yrs to ensure same random number sequence for 
		  observation errors.
 	PSUEDOCODE:
 		1)  calculate selectivities to be used in the simulations.
 		2)  calculate mortality rates (M), F is conditioned on catch.
 		3)  generate random numbers for observation & process errors.
 		4)  calculate survivorship and stock-recruitment pars based on average M & fec.
 		5)  initialize state variables.
 		6)  population dynamics with F conditioned on catch.
 		7)  compute catch-at-age samples.
 		8)  compute total catch.
 		9)  compute the relative abundance indices.
 		10) rewrite empitical weight at age matrix.
 		11) write simulated data to file.
  	TODO list:
	[?] - March 9, 2013.  Fix simulation model to generate consistent data when 
		  doing retrospective analyses on simulated datasets.
	[ ] - TODO: Stock-Recruitment model.
	[ ] - TODO: switch statement for catch-type to get Fishing mortality rate.
  	[ ] 
  	*/
}

void model_parameters::simulationModel(const long& seed)
{
  {
	cout<<global_parfile<<endl;
	bool pinfile = 0;
	cout<<"___________________________________________________\n"<<endl;
	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
	cout<<"___________________________________________________"<<endl;
	//if(norm(log_rec_devs)!=0)
	if( global_parfile )
	{
		cout<<"\tUsing pin file for simulation"<<endl;
		pinfile = 1;
	}
	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
	cout<<"\tNumber of retrospective years: "<<retro_yrs<<endl;
	cout<<"___________________________________________________\n"<<endl;
	int ii,ig,ih;
	// |---------------------------------------------------------------------------------|
	// | 1) SELECTIVITY
	// |---------------------------------------------------------------------------------|
	// |
    calcSelectivities(isel_type);
    // |---------------------------------------------------------------------------------|
    // | 2) MORTALITY
    // |---------------------------------------------------------------------------------|
    // | - NOTE only natural mortality is computed at this time.
    // | [ ] - add simulated random-walk in natural mortality rate here.
    // |
    calcTotalMortality();
    F.initialize();
    Z.initialize();
    S.initialize();
    // |---------------------------------------------------------------------------------|
    // | 3) GENERATE RANDOM NUMBERS
    // |---------------------------------------------------------------------------------|
    // | - epsilon -> Observation errors
    // | - rec_dev -> Process errors
    // | - init_rec_dev
    // | [ ] - add other required random numbers if necessary.
    // |
	random_number_generator rng(seed);
	dmatrix      epsilon(1,nItNobs,1,n_it_nobs);
	dmatrix      rec_dev(1,n_ag,syr,nyr+retro_yrs);
	dmatrix init_rec_dev(1,n_ag,sage+1,nage);
	dvector      eta(1,nCtNobs);
	epsilon.fill_randn(rng);
	rec_dev.fill_randn(rng);
	init_rec_dev.fill_randn(rng);
	eta.fill_randn(rng);
    // | Scale survey observation errors
    double std;
    for(k=1;k<=nItNobs;k++)
    {
    	for(i=1;i<=n_it_nobs(k);i++)
    	{
    		std = 1.0e3;
    		if( it_wt(k,i)>0 )
    		{
    			std = value(sig(it_grp(k,i))/it_wt(k,i));
    		}
    		epsilon(k,i) = epsilon(k,i)*std - 0.5*std*std;
    	}
    }
    // | Scale process errors
    for(ih=1;ih<=n_ag;ih++)
    {
		std              = value(tau(1));
		rec_dev(ih)      = rec_dev(ih) * std - 0.5*std*std;
		init_rec_dev(ih) = init_rec_dev(ih)*std - 0.5*std*std;
    }
    // | Scale total catch errors
    std = d_iscamCntrl(4);
    for(ii=1;ii<=nCtNobs;ii++)
    {
    	eta(ii) = eta(ii)* std  - 0.5*std*std;
    }
    // |---------------------------------------------------------------------------------|
    // | 4) SURVIVORSHIP & STOCK-RECRUITMENT PARAMETERS BASED ON AVERAGE M & FECUNDITY
    // |---------------------------------------------------------------------------------|
    // | -> Loop over each group/stock and compute survivorship, phib, so and beta.
    // | - fa is the average mature weight-at-age
    // |
    double phib;
    dvector ma(sage,nage);
    dvector fa(sage,nage);
    dvector lx(sage,nage);
    dvector lw(sage,nage);
    for(g=1;g<=ngroup;g++)
    {
    	lx.initialize();
		lw.initialize();
		lx(sage) = 1.0;
		lw(sage) = 1.0;
		phib = 0;
		for(f=1;f<=narea;f++)
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);
				for(j=sage;j<=nage;j++)
				{
					ma(j) = mean( trans(value(M(ig)))(j)  );
					fa(j) = mean( trans(d3_wt_mat(ig))(j) );
					if(j>sage)
					{
						lx(j) = lx(j-1) * exp(-ma(j-1));
					}
					lw(j) = lx(j) * exp( -ma(j)*d_iscamCntrl(13) );
				}
				lx(nage) /= 1.0 - exp(-ma(nage));
				lw(nage) /= 1.0 - exp(-ma(nage));
				phib += 1./(narea*nsex) * lw * fa;
			}
		}
		so(g)  = kappa(g)/phib;
		sbo(g) = ro(g) * phib;
		switch(int(d_iscamCntrl(2)))
		{
			case 1:
				beta(g) = (kappa(g)-1.0)/sbo(g);
			break;
			case 2:
				beta(g) = log(kappa(g))/sbo(g);
			break;
		}
    }
   // |---------------------------------------------------------------------------------|
   // | 5) INITIALIZE STATE VARIABLES
   // |---------------------------------------------------------------------------------|
   // |
	N.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		dvector tr(sage,nage);
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);
		lx.initialize();
		lx(sage) = 1.0;
		for(j=sage;j<nage;j++)
		{
			lx(j+1) = lx(j) * exp(-value(M(ig)(syr)(j)));
		}
		lx(nage) /= 1.0 - exp(-value(M(ig)(syr)(nage)));
		if( d_iscamCntrl(5) )
		{
			tr = log( value(ro(g)) ) + log(lx);
		}
		else if( !d_iscamCntrl(5) )
		{
			tr(sage)        = value(log_avgrec(ih)+rec_dev(ih)(syr));;
			tr(sage+1,nage) = value(log_recinit(ih)+init_rec_dev(ih));
			tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
		}
		N(ig)(syr)(sage,nage) = 1./nsex * exp(tr);
		log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);
		for(i=syr+1;i<=nyr;i++)
		{
			log_rt(ih)(i) = (log_avgrec(ih)+rec_dev(ih)(i));
			N(ig)(i,sage) = 1./nsex * exp( log_rt(ih)(i) );
		}
		N(ig)(nyr+1,sage) = 1./nsex * exp(log_avgrec(ih));
	}
	// |---------------------------------------------------------------------------------|
	// | 6) POPULATION DYNAMICS WITH F CONDITIONED ON OBSERVED CATCH
	// |---------------------------------------------------------------------------------|
	// | - va  -> matrix of fisheries selectivity coefficients.
	// | - [ ] TODO: switch statement for catch-type to get Fishing mortality rate.
	// | - [ ] TODO: Stock-Recruitment model (must loop over area sex for each group).
	// | - bug! ft is a global variable used in calcCatchAtAge and calcTotalCatch. FIXED
	// | - [ ] TODO: fishing mortality rate with sex ratio of catch is unknown and assume
	// |             the same F value for both male and females.
	dmatrix va(1,ngear,sage,nage);
	dmatrix zt(syr,nyr,sage,nage);
	dmatrix st(syr,nyr,sage,nage);
	BaranovCatchEquation cBaranov;
	for(ig=1;ig<=n_ags;ig++)
	{
		dmatrix tmp_ft(syr,nyr,1,ngear);
		for(i=syr;i<=nyr;i++)
		{
			dvector ba = elem_prod(value(N(ig)(i)),d3_wt_avg(ig)(i));
			dvector ct = d3_Ct(ig)(i);
			// | Selectivity modifications if necessary
			for(k=1;k<=ngear;k++)
			{
				va(k) = exp(value(log_sel(k)(ig)(i)));
				if( d_iscamCntrl(15) == 1 && dAllocation(k) > 0 )
				{
					va(k)             = ifdSelex(va(k),ba,0.25);
					log_sel(k)(ig)(i) = log(va(k));
					// dlog_sel(k)(i) = log(va(k));
				}
			}
			// | [ ] TODO switch statement for catch_type to determine F.
			tmp_ft(i) = cBaranov.getFishingMortality(ct, value(M(ig)(i)), va, value(N(ig)(i)),d3_wt_avg(ig)(i));
			zt(i) = value(M(ig)(i));
			for(k=1;k<=ngear;k++)
			{
				ft(ig)(k)(i) = tmp_ft(i,k);
				F(ig)(i) += tmp_ft(i,k) * va(k);
				zt(i) += tmp_ft(i,k) * va(k);
			}
			st(i) = exp(-zt(i));
			Z(ig)(i) = M(ig)(i) + F(ig)(i);
			S(ig)(i) = exp(-Z(ig)(i));
			// | [ ] TODO: Stock-Recruitment model
			// | [ ] TODO: Add autocorrelation.
			if( !pinfile ) COUT("Add stock recruitment Model")
			/* 
			sbt(ig)(i) = (elem_prod(N(ig)(i),exp(-zt(i)*d_iscamCntrl(13)))*d3_wt_mat(ig)(i));
			if(i>=syr+sage-1 && !pinfile)
			{
				double rt,et;
				et = value(sbt(ig)(i-sage+1));
				if(d_iscamCntrl(2)==1)
				{
					rt = value(so(g)*et/(1.+beta(g)*et));
				}
				else
				{
					rt = value(so(g)*et*exp(-beta(g)*et));
				}
				N(ig)(i)(sage) = rt * exp(rec_dev(ih)(i)-0.5*tau*tau);
			}
			*/
			// | Update state variables
			N(ig)(i+1)(sage+1,nage) =++ elem_prod(N(ig)(i)(sage,nage-1),st(i)(sage,nage-1));
			N(ig)(i+1)(nage) += N(ig)(i)(nage)*st(i)(nage);
		}
	}
	// |---------------------------------------------------------------------------------|
	// | 7) CATCH-AT-AGE
	// |---------------------------------------------------------------------------------|
	// | - A is the matrix of observed catch-age data.
	// | - A_hat is the predicted matrix from which to draw samples.
	// |
	int k, kk, aa, AA;
	double age_tau = value(sig(1));
	calcComposition();
	for(kk=1;kk<=nAgears;kk++)
	{
		aa = n_A_sage(kk);
		AA = n_A_nage(kk);
		dvector pa(aa,AA);
		for(ii=1;ii<=n_A_nobs(kk);ii++)
		{
			pa = value(A_hat(kk)(ii));
			d3_A(kk)(ii)(aa,AA)=rmvlogistic(pa,age_tau,i+seed);
		}
	}
	// |---------------------------------------------------------------------------------|
	// | 8) TOTAL CATCH
	// |---------------------------------------------------------------------------------|
	// | - dCatchData is the matrix of observations
	// | - need to over-write the d3_Ct with the new errors.
	calcTotalCatch();
	d3_Ct.initialize();
	for(ii=1;ii<=nCtNobs;ii++)
	{
		dCatchData(ii,7) = value(ct(ii)) * exp(eta(ii));
		i = dCatchData(ii)(1);
		k = dCatchData(ii)(2);
		f = dCatchData(ii)(3);
		g = dCatchData(ii)(4);
		h = dCatchData(ii)(5);
		if( h==0 )
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);
				d3_Ct(ig)(i)(k) = 1./nsex*dCatchData(ii)(7);
			}
		}
		else
		{
			ig = pntr_ags(f,g,h);
			d3_Ct(ig)(i)(k) = dCatchData(ii)(7);
		} 
	}
	// cout<<d3_Ct(1)<<endl;
	// |---------------------------------------------------------------------------------|
	// | 9) RELATIVE ABUNDANCE INDICES
	// |---------------------------------------------------------------------------------|
	// | - d3_survey_data is the matrix of input data.
	// |
	calcSurveyObservations();
	for(kk=1;kk<=nItNobs;kk++)
	{
		for(ii=1;ii<=n_it_nobs(kk);ii++)
		{
			d3_survey_data(kk)(ii)(2) *= exp(epsilon(kk)(ii));	
		}
	}
	// |---------------------------------------------------------------------------------|
	// | 10) Empirical weight at age
	// |---------------------------------------------------------------------------------|
	// | 
		for(int ii=1; ii<=nWtTab; ii++)
		{
			if(nWtNobs(ii) > 0 )
			{
				d3_inp_wt_avg(ii,1,sage-5) = d3_inp_wt_avg(ii,1,sage-5)*projwt(ii);
			}
		}
		cout<< d3_inp_wt_avg(1)(1)(sage-5,nage) <<endl;
	// |---------------------------------------------------------------------------------|
	// | 11) WRITE SIMULATED DATA TO FILE
	// |---------------------------------------------------------------------------------|
	// |
	writeSimulatedDataFile();
  }
  	/**
  	Purpose:  This function writes a simulated data file based on the simulation
  		  model output when the user specifies the -sim option.  This is only
  	      necessary if the user wishes to perform a retrospecrtive analysis on
  	      simulated data. 
  	Author: Steven Martell
  	Arguments:
  	\param	seed -> the random number seed that is concatenated into the file name.
  	NOTES:
  	TODO list:
  	[ ] 
  	*/
}

void model_parameters::writeSimulatedDataFile(void)
{
  {
  	adstring sim_datafile_name = "Simulated_Data_"+str(rseed)+".dat";
  	ofstream dfs(sim_datafile_name);
  	dfs<<"#Model dimensions"<<endl;
  	dfs<< narea 		<<endl;
  	dfs<< ngroup		<<endl;
  	dfs<< nsex			<<endl;
  	dfs<< syr   		<<endl;
  	dfs<< nyr   		<<endl;
  	dfs<< sage  		<<endl;
  	dfs<< nage  		<<endl;
  	dfs<< ngear 		<<endl;
  	dfs<<"#Allocation"	<<endl;
  	dfs<< dAllocation 	<<endl;
  	dfs<<"#Age-schedule and population parameters"<<endl;
  	dfs<< d_linf  			<<endl;
  	dfs<< d_vonbk  			<<endl;
  	dfs<< d_to  			<<endl;
  	dfs<< d_a  				<<endl;
  	dfs<< d_b  				<<endl;
  	dfs<< d_ah  			<<endl;
  	dfs<< d_gh  			<<endl;
  	dfs<< n_MAT				<<endl;
	dfs<< d_maturityVector <<endl;
  	dfs<<"#Observed catch data"<<endl;
  	dfs<< nCtNobs 		<<endl;
  	dfs<< dCatchData    <<endl;
  	dfs<<"#Abundance indices"	<<endl;
  	dfs<< nItNobs 					<<endl;
  	dfs<< n_it_nobs 				<<endl;
  	dfs<< n_survey_type 			<<endl;
  	dfs<< d3_survey_data 			<<endl;
  	dfs<<"#Age composition"		<<endl;
  	dfs<< nAgears				<<endl;
  	dfs<< n_A_nobs				<<endl;
  	dfs<< n_A_sage				<<endl;
  	dfs<< n_A_nage				<<endl;
  	dfs<< inp_nscaler 			<<endl;
  	dfs<< d3_A					<<endl;
  	dfs<<"#Empirical weight-at-age data"	<<endl;
  	dfs<< nWtTab 				<<endl;
  	dfs<< nWtNobs				<<endl;
	dfs<< d3_inp_wt_avg			<<endl; // not sure if this shoud be d3_inp_wt_avg, and how this would affect simDatfile 
	dfs<<"#EOF"	<<endl;
	dfs<< 999	<<endl;
	// | END OF WRITING SIMULATED DATAFILE.
  }
  	/**
  	Purpose:  This function returns a modified selectivity vector (va) based on
  			  the assumption that age-based selectivity will operate on the principle
  	          of ideal free distribution.
  	Author: Steven Martell
  	Arguments:
  	\param	va -> age-specific vulnerability
  	\param	ba -> age-specific biomass (relative abundance is fine)
  	NOTES:
  	TODO list:
  	[ ] 
  	*/
}

dvector model_parameters::ifdSelex(const dvector& va, const dvector& ba, const double& mpow)
{
  {
  	dvector pa(sage,nage);
  	pa = (elem_prod(va,pow(ba,mpow)));
  	pa = pa/sum(pa);
  	pa = exp( log(pa) - log(mean(pa)) );
  	return (pa);
  }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	if(verbose)cout<<"Start of Report Section..."<<endl;
	report<<DataFile<<endl;
	report<<ControlFile<<endl;
	report<<ProjectFileControl<<endl;
	REPORT(objfun);
	REPORT(nlvec);
	REPORT(ro);
	dvector rbar=value(exp(log_avgrec));
	REPORT(rbar);
	dvector rinit=value(exp(log_recinit));
	REPORT(rinit);
	REPORT(sbo);
	REPORT(kappa);
	dvector steepness=value(theta(2));
	REPORT(steepness);
	REPORT(m);
	// double tau = value(sqrt(1.-rho)*varphi);
	// double sig = value(sqrt(rho)*varphi);
	report<<"rho"<<endl<<theta(6)<<endl;
	report<<"vartheta"<<endl<<theta(7)<<endl;
	REPORT(varphi);
	REPORT(tau);
	REPORT(sig);
	REPORT(age_tau2);
	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS & AGE-SCHEDULE INFORMATION ON GROWTH AND MATURITY
	// |---------------------------------------------------------------------------------|
	// |
	REPORT(narea);
	REPORT(ngroup);
	REPORT(nsex);
	REPORT(syr);
	REPORT(nyr);
	REPORT(sage);
	REPORT(nage);
	REPORT(ngear);
	ivector yr(syr,nyr);
	ivector yrs(syr,nyr+1);
	yr.fill_seqadd(syr,1); 
	yrs.fill_seqadd(syr,1); 
	REPORT(yr);
	REPORT(yrs);
	// REPORT(iyr);  //DEPRECATE, old survey years
	REPORT(age);
	REPORT(la);
	REPORT(wa);
	REPORT(ma);
	// |---------------------------------------------------------------------------------|
	// | OBSERVED AND PREDICTED DATA AND RESIDUALS
	// |---------------------------------------------------------------------------------|
	// | - Catch data
	// | - Survey data
	// | - Age composition data
	// | - Empirical weight-at-age data
	REPORT(dCatchData);
	REPORT(ct);
	REPORT(eta);
	REPORT(q);
	REPORT(qt);
	REPORT(d3_survey_data);
	REPORT(it_hat);
	REPORT(it_wt);
	REPORT(epsilon);
	//RF changed this to 1
	if(n_A_nobs(1) > 0)
	{
		REPORT(n_A_sage);
		REPORT(n_A_nage);
		REPORT(d3_A);
		REPORT(A_hat);
		REPORT(A_nu);
		/// The following is a total hack job to get the effective sample size
		/// for the multinomial distributions.
		// FIXED the retrospective bug here near line 4507 (if iyr<=nyr)
		report<<"Neff"<<endl;
		dvector nscaler(1,nAgears);
		nscaler.initialize();
   	 int naa;
   	 int iyr;
		for(k = 1; k<=nAgears; k++)
		{
			if( int(nCompLikelihood(k)) )
			{
				naa = 0;
				//retrospective counter
				for(i=1;i<=n_A_nobs(k);i++)
				{
					iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
					if(iyr<=nyr) naa++; else continue;
				}
				dmatrix     O = trans(trans(d3_A_obs(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
				dvar_matrix P = trans(trans(A_hat(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
				for(j = 1; j<= naa; j++)
				{
					double effectiveN = neff(O(j)/sum(O(j)),P(j));
					report<<sum(O(j))<<"\t"<<effectiveN<<endl;
					nscaler(k) += effectiveN;
				}
				nscaler(k) /= naa;
			}
		}
		REPORT(nscaler);
	}
	// d3_wt_avg(1,n_ags,syr,nyr+1,sage,nage);
	adstring tt = "\t";
	//REPORT(xxinp_wt_avg);
	REPORT(dWt_bar);
	// REPORT(d3_wt_avg);
	REPORT(d3_wt_mat);
	REPORT(d3_wt_dev);
	report<<"d3_wt_avg"<<endl;
	for(int ig=1;ig<=n_ags;ig++)
	{
		f = n_area(ig);
		g = n_group(ig);
		h = n_sex(ig);
		for(i=syr;i<=nyr;i++)
		{
			//year area stock sex |age columns (sage, nage) of weight at age data |
			report<<i<<tt;
			report<<f<<tt;
			report<<g<<tt;
			report<<h<<tt;
			report<<d3_wt_avg(ig)(i)<<endl;
		}
	}
	// |---------------------------------------------------------------------------------|
	// | SELECTIVITIES (4darray)
	// |---------------------------------------------------------------------------------|
	// |
	report<<"log_sel"<<endl;
	for(k=1;k<=ngear;k++)
	{
		for(int ig=1;ig<=n_ags;ig++)
		{
			for(i=syr;i<=nyr;i++)
			{
				report<<k<<"\t"<<ig<<"\t"<<i<<"\t"<<log_sel(k)(ig)(i)<<endl;	
			}
		}
	}
	// |---------------------------------------------------------------------------------|
	// | MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	// REPORT(ft);
	report<<"ft"<<endl;
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		report<<ft(ig)<<endl;
	}
	REPORT(M);
	REPORT(F);
	REPORT(Z);
	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT
	// |---------------------------------------------------------------------------------|
	// |
	int rectype=int(d_iscamCntrl(2));
	REPORT(rectype);
	REPORT(so);
	REPORT(beta);
	REPORT(sbt);
	REPORT(bt);
	REPORT(rt);
	REPORT(delta);
	dmatrix rep_rt = value( exp(trans(trans(log_rt).sub(syr,nyr))) );
	for(int ig = 1; ig <= n_ag; ig++ )
	{
		rep_rt(ig)(syr) = value( exp( log_rt(ig)(syr-nage+sage) ) );
	}
	REPORT(rep_rt);
	// |---------------------------------------------------------------------------------|
	// | ABUNDANCE IN NUMBERS 
	// |---------------------------------------------------------------------------------|
	// |
	REPORT(N);
	//START_RF_ADD
	// |---------------------------------------------------------------------------------|
	// | ANNUAL MEAN WEIGHT DATA 
	// |---------------------------------------------------------------------------------|
	// |
	REPORT(obs_annual_mean_weight);
	REPORT(annual_mean_weight);
	//END_RF_ADD
	// |---------------------------------------------------------------------------------|
	// | MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	if( last_phase() )
	{
		cout<<"Calculating MSY-based reference points"<<endl;
		calcReferencePoints();
		cout<<"Finished calcReferencePoints"<<endl;
		//exit(1);
		REPORT(bo);
		REPORT(fmsy);
		REPORT(msy);
		REPORT(bmsy);
		// REPORT(Umsy);
	}
	cout<<"You got to the end of the report section"<<endl;
	// |---------------------------------------------------------------------------------|
	// | OUTPUT FOR OPERATING MODEL
	// |---------------------------------------------------------------------------------|
	// | Move to final section?
	if( last_phase() )
	{
		ofstream ofs("iSCAM.res");
		ofs<<"# Bo\n"<<bo<<endl;
		ofs<<"# Fmsy\n"<<fmsy<<endl;
		ofs<<"# MSY\n"<<msy<<endl;
		ofs<<"# Bmsy\n"<<bmsy<<endl;
		ofs<<"# Sbt\n";
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<sbt(g)(nyr+1)<<"\t";
		}
		ofs<<endl;
		// projected biomass
		// The total biomass for each stock
		ofs<<"# Total biomass\n";
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<bt(g)(nyr+1)<<"\t";
		}
		ofs<<endl;
		ofs<<"# Numbers-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<N(ig)(nyr+1)<<endl;
		}
		ofs<<"# Weight-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<d3_wt_avg(ig)(nyr+1)<<endl;
		}		
		ofs<<"# Natural mortality-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<M(ig)(nyr)<<endl;
		}		
		// 4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
		ofs<<"# log_selectivity\n";
		for(int k = 1; k <= ngear; k++ )	
		{
			for(int ig = 1; ig <= n_ags; ig++ )
			{
				ofs<<log_sel(k)(ig)(nyr)<<endl;
			}
		}
	}
	if(verbose)cout<<"END of Report Section..."<<endl;
	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	//#if defined __APPLE__ || defined __linux
	//if( retro_yrs && last_phase() )
	//{
	//	//adstring rep="iscam.ret"+str(retro_yrs);
	//	//rename("iscam.rep",rep);
	//	adstring copyrep = "cp iscam.rep iscam.ret"+str(retro_yrs);
	//	system(copyrep);
	//}
	//// cout<<"Ya hoo"<<endl;
	////exit(1);
	//#endif
	//#endif
   // REPORT_SECTION END
}

void model_parameters::generate_new_files(void)
{
  ofstream rd("RUN.dat");
  rd<<NewFileName + ".dat"<<endl;
  rd<<NewFileName + ".ctl"<<endl;
  rd<<NewFileName + ".pfc"<<endl;
  exit(1);
  #if defined __APPLE__ || defined __linux
    adstring bscmddat = "cp ../lib/iscam.dat" + NewFileName +".dat";
    system(bscmddat);
    adstring bscmdctl = "cp ../lib/ iscam.ctl" + NewFileName +".ctl";
    system(bscmdctl);
    adstring bscmdpfc = "cp ../lib/ iscam.PFC" + NewFileName +".pfc";
    system(bscmdpfc);	
  #endif
  #if defined _WIN32 || defined _WIN64
    adstring bscmddat = "copy ../lib/iscam.dat" + NewFileName +".dat";
    system(bscmddat);
    adstring bscmdctl = "copy ../lib/ iscam.ctl" + NewFileName +".ctl";
    system(bscmdctl);
    adstring bscmdpfc = "copy ../lib/ iscam.PFC" + NewFileName +".pfc";
    system(bscmdpfc);	
  #endif
}

void model_parameters::mcmc_output(void)
{
  if(nf==1){
    // Open the files and write the headers
    ofstream ofs("iscam_mcmc.csv");
    // The structure for these objects can be found at roughly lines 924 and 1409.
    // they are set up as vector_vectors to increase dimensionality
    // for the split sex case and also for areas and groups
    // The format for the mcmc output headers will be one of:
    // parametername_gr[0-9]+  - for unique group only
    // parametername_gs[0-9]+  - for unique group and sex
    // paramatername_ag[0-9]+  - for unique area and gear
    for(int group=1;group<=ngroup;group++){
      ofs<<"ro_gr"<<group;
    }
    for(int group=1;group<=ngroup;group++){
      ofs<<","<<"h_gr"<<group;
    }
    for(int gs=1;gs<=n_gs;gs++){
      ofs<<","<<"m_gs"<<gs;
    }
    for(int ag=1;ag<=n_ag;ag++){
      ofs<<","<<"rbar_ag"<<ag;
    }
    for(int ag=1;ag<=n_ag;ag++){
      ofs<<","<<"rinit_ag"<<ag;
    }
    for(int group=1;group<=ngroup;group++){
      ofs<<","<<"rho_gr"<<group;
    }
    for(int group=1;group<=ngroup;group++){
      ofs<<","<<"vartheta_gr"<<group;
    }
    ofs<<","<<"bo";
    ofs<<","<<"bmsy";
    for(int fleet=1;fleet<=nfleet;fleet++){
      ofs<<","<<"msy"<<fleet;
    }
    for(int fleet=1;fleet<=nfleet;fleet++){
      ofs<<","<<"fmsy"<<fleet;
    }
    for(int i=1;i<=nItNobs;i++){
      ofs<<","<<"q"<<i;
    }
    for(int group=1;group<=ngroup;group++){
      ofs<<","<<"SSB"<<group;
    }
    ofs<<","<<"f";
    ofs<<endl;
    ofstream of1("iscam_sbt_mcmc.csv");
    for(int group=1;group<=ngroup;group++){
      for(int yr=syr;yr<=nyr+1;yr++){
        if(yr == syr){
          of1<<"sbt"<<group<<"_"<<yr;
        }else{
          of1<<",sbt"<<group<<"_"<<yr;
        }
      }
    }
    of1<<endl;
    ofstream of2("iscam_rt_mcmc.csv");
    for(int group=1;group<=ngroup;group++){
      for(int yr=syr+sage;yr<=nyr;yr++){
        if(yr == syr+sage){
          of2<<"rt"<<group<<"_"<<yr;
        }else{
          of2<<",rt"<<group<<"_"<<yr;
        }
      }
    }
    of2<<endl;
    ofstream of3("iscam_ft_mcmc.csv");
    int iter = 1;
    for(int ag=1;ag<=n_ags;ag++){
      for(int gear=1;gear<=ngear;gear++){
        for(int yr=syr;yr<=nyr;yr++){
          if(iter == 1){
            of3<<"ft"<<ag<<"_gear"<<gear<<"_"<<yr;
          }else{
            of3<<",ft"<<ag<<"_gear"<<gear<<"_"<<yr;
          }
          iter++;
        }
      }
    }
    of3<<endl;
  }
  // Leading parameters & reference points
  calcReferencePoints();
  // Append the values to the files
  ofstream ofs("iscam_mcmc.csv",ios::app);
  for(int group=1;group<=ngroup;group++){
    ofs<<exp(theta(1)(group));
  }
  for(int group=1;group<=ngroup;group++){
    ofs<<","<<theta(2)(group);
  }
  for(int gs=1;gs<=n_gs;gs++){
    ofs<<","<<exp(theta(3)(gs));
  }
  for(int ag=1;ag<=n_ag;ag++){
    ofs<<","<<exp(theta(4)(ag));
  }
  for(int ag=1;ag<=n_ag;ag++){
    ofs<<","<<exp(theta(5)(ag));
  }
  for(int group=1;group<=ngroup;group++){
    ofs<<","<<theta(6)(group);
  }
  for(int group=1;group<=ngroup;group++){
    ofs<<","<<theta(7)(group);
  }
  ofs<<","<<bo;
  ofs<<","<<bmsy;
  for(int fleet=1;fleet<=nfleet;fleet++){
    ofs<<","<<msy(fleet);
  }
  for(int fleet=1;fleet<=nfleet;fleet++){
    ofs<<","<<fmsy(fleet);
  }
  for(int it=1;it<=nItNobs;it++){
    ofs<<","<<q(it);
  }
  for(int group=1;group<=ngroup;group++){
    ofs<<","<<sbt(group)(nyr);
  }
  ofs<<","<<objfun;
  ofs<<endl;
  // output spawning stock biomass
  ofstream of1("iscam_sbt_mcmc.csv",ios::app);
  for(int group=1;group<=ngroup;group++){
    for(int yr=syr;yr<=nyr+1;yr++){
      if(yr == syr){
        of1<<sbt(group)(yr);
      }else{
        of1<<","<<sbt(group)(yr);
      }
    }
  }
  of1<<endl;
  // output age-1 recruits
  ofstream of2("iscam_rt_mcmc.csv",ios::app);
  for(int group=1;group<=ngroup;group++){
    for(int yr=syr+sage;yr<=nyr;yr++){
      if(yr == syr+sage){
        of2<<rt(group)(yr);
      }else{
        of2<<","<<rt(group)(yr);
      }
    }
  }
  of2<<endl;
  // output fishing mortality
  ofstream of3("iscam_ft_mcmc.csv",ios::app);
  int iter = 1;
  for(int ag=1;ag<=n_ags;ag++){
    for(int gear=1;gear<=ngear;gear++){
      for(int yr=syr;yr<=nyr;yr++){
        if(iter == 1){
          of3<<ft(ag)(gear)(yr);
        }else{
          of3<<","<<ft(ag)(gear)(yr);
        }
        iter++;
      }
    }
  }
  of3<<endl;
  ofs.flush();
  of1.flush();
  of2.flush();
  of3.flush();
}

void model_parameters::runMSE()
{
	cout<<"Start of runMSE"<<endl;
	// STRUCT FOR MODEL VARIABLES
	ModelVariables s_mv;
	s_mv.log_ro    = value( theta(1) );
	s_mv.steepness = value( theta(2) );
	s_mv.m         = value( theta(3) );
	s_mv.log_rbar  = value( theta(4) );
	s_mv.log_rinit = value( theta(5) );
	s_mv.rho       = value( theta(6) );
	s_mv.varphi    = value( theta(7) );
	// Selectivity parameters
	d3_array log_sel_par(1,ngear,1,jsel_npar,1,isel_npar);
	d4_array d4_log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
	for(int k = 1; k <= ngear; k++ )
	{
		log_sel_par(k) = value(sel_par(k));
		d4_log_sel(k)  = value(log_sel(k));
	}
	s_mv.d3_log_sel_par = &log_sel_par;
	s_mv.d4_logSel      = &d4_log_sel;
	d3_array d3_M(1,n_ags,syr,nyr,sage,nage);
	d3_array d3_F(1,n_ags,syr,nyr,sage,nage);
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		d3_M(ig) = value(M(ig));
		d3_F(ig) = value(F(ig));
	}
	s_mv.d3_M = &d3_M;
	s_mv.d3_F = &d3_F;
	s_mv.log_rec_devs = value(log_rec_devs);
	s_mv.init_log_rec_devs = value(init_log_rec_devs);
	s_mv.q = value(q);
	s_mv.sbt = value(sbt);
	d3_array tmp_ft=value(ft);
	s_mv.d3_ft = &tmp_ft;
	s_mv.sbo = value(sbo);
	s_mv.so = value(so);
	// |-----------------------------------|
	// | Instantiate Operating Model Class |
	// |-----------------------------------|
	OperatingModel om(s_mv,argc,argv);
	om.runScenario(rseed);
	COUT("DONE");
}

void model_parameters::final_calcs()
{
	cout<<"Here I am "<<endl;
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"--Number of function evaluations: "<<nf<<endl;
	cout<<"--Results are saved with the base name:\n"<<"\t"<<BaseFileName<<endl;
	cout<<"*******************************************"<<endl;
	if(mseFlag) runMSE();
	cout<<"End of class testing"<<endl;
	//Make copies of the report file using the ReportFileName
	//to ensure the results are saved to the same directory 
	//that the data file is in. This should probably go in the 
	//FINAL_SECTION
	//CHANGED only copy over the mcmc files if in mceval_phase()
	#if defined __APPLE__ || defined __linux
	if(last_phase() && !retro_yrs)
	{
		adstring bscmd = "cp iscam.rep " +ReportFileName;
		system(bscmd);
		bscmd = "cp iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		bscmd = "cp iscam.std " + BaseFileName + ".std";
		system(bscmd);
		bscmd = "cp iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
		//if( SimFlag )
		//{
		//	bscmd = "cp iscam.sim " + BaseFileName + ".sim";
		//	system(bscmd);
		//}
			ofstream mcofs(ReportFileName,ios::app);
			mcofs<<"ENpar\n"<<dicNoPar<<endl;
			mcofs<<"DIC\n"<<dicValue<<endl;
			mcofs.close();
			cout<<"Copied MCMC Files"<<endl;
		}
	if( last_phase() && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "cp iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}
	#endif
	#if defined _WIN32 || defined _WIN64
	if(last_phase() && !retro_yrs)
	{
		adstring bscmd = "copy iscam.rep " +ReportFileName;
		system(bscmd);
		bscmd = "copy iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		bscmd = "copy iscam.std " + BaseFileName + ".std";
		system(bscmd);
		bscmd = "copy iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
	}
	if( last_phase() && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "copy iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}
	#endif
	if(mseFlag) runMSE();
	cout<<"End of class testing"<<endl;
	//exit(1);
	//  Print run time statistics to the screen.
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"--Number of function evaluations: "<<nf<<endl;
	cout<<"--Results are saved with the base name:\n"<<"\t"<<BaseFileName<<endl;
	cout<<"*******************************************"<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
