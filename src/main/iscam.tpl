/// iSCAM = integrated Statistical Catch Age Model

DATA_SECTION
	// |---------------------------------------------------------------------------------|
	// | STRINGS FOR INPUT FILES                                                         |
	// |---------------------------------------------------------------------------------|
	/// | DataFile.dat           : data to condition the assessment model on
	init_adstring DataFile;      ///< String for the input datafile name.
	/// | ControlFile.ctl        : controls for phases, selectivity options
	init_adstring ControlFile;   ///< String for the control file.
	/// | ProjectFileControl.pfc : used for stock projections under TAC
	init_adstring ProjectFileControl;  ///< String for the projection file.

	init_adstring ProcedureControlFile;

	init_adstring ScenarioControlFile;
	/// | BaseFileName           : file prefix used for all iSCAM model output
	!! BaseFileName = stripExtension(ControlFile);  ///< BaseName given by the control file
	/// | ReportFileName         : file name to copy report file to.
	!! ReportFileName = BaseFileName + adstring(".rep");
	!! LOG<<"Base part of file names:\n"<<BaseFileName<<'\n';

  // |---------------------------------------------------------------------------------|
	// | READ IN PROJECTION FILE CONTROLS
	// |---------------------------------------------------------------------------------|
	// | n_tac    : length of catch vector for decision table catch stream.
	// | tac      : vector of total catch values to be used in the decision table.
	// | pf_cntrl : vector of controls for the projections.
	// | Documentation for projection control file pf_cntrl
	// | 1) start year for m_bar calculation
	// | 2)   end year for m_bar calculation
	// | 3) start year for average fecundity/weight-at-age
	// | 4)   end year for average fecundity/weight-at-age
	// | 5) start year for recruitment period (not implemented yet)
	// | 6)   end year for recruitment period (not implemented yet)
	// |

	!! ad_comm::change_datafile_name(ProjectFileControl);
	/// | Number of catch options to explore in the decision table.
	init_int n_tac; ///< Number of catch options to explore in the decision table.
  !! LOG<<"Number of TAC options to explore in decision table:\n";
  !! LOG<<n_tac<<'\n';
	/// | Vector of catch options.
	init_vector tac(1,n_tac);
	init_int n_pfcntrl;
	init_vector pf_cntrl(1,n_pfcntrl);


	//init_vector mse_cntrl(1,1);

	init_int eof_pf;

	LOC_CALCS
		if(eof_pf!=-999)
		{
			LOG<<"Error reading projection file.\n";
			LOG<<"Last integer read is "<<eof_pf<<'\n';
			LOG<<"The file should end with -999.\n Aborting!\n";
			ad_exit(1);
		}
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | COMMAND LINE ARGUMENTS FOR SIMULATION & RETROSPECTIVE ANALYSIS
	// |---------------------------------------------------------------------------------|
	// | SimFlag    : if user specifies -sim, then turn SimFlag on.
	// | retro_yrs  : number of terminal years to remove.

	int SimFlag;  ///< Flag for simulation mode
	int mseFlag;  ///< Flag for management strategy evaluation mode
	int rseed;    ///< Random number seed for simulated data.
	int retro_yrs;///< Number of years to look back from terminal year.
	int testMSY;

	int delaydiff; ///Flag for delay difference model 

	LOC_CALCS
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

		// command line option for retrospective analysis. "-retro retro_yrs"
		retro_yrs=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retro_yrs = atoi(ad_comm::argv[on+1]);
			LOG<<" ____________________________________________________ \n";
			LOG<<"| Implementing Retrospective analysis                |\n";
			LOG<<"|____________________________________________________|\n";
			LOG<<"| Number of retrospective years = "<<retro_yrs<<'\n';
		}

		// Management strategy evaluation.
		mseFlag = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt))>-1)
		{
			mseFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
			LOG<<" _________________________________________________ \n";
			LOG<<"|Implementing Management Strategy Evaluation      |\n";
			LOG<<"|_________________________________________________|\n";
		}

		// Test MSY
		testMSY = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-msy",opt))>-1){
			LOG<<"Testing MSY calculations with Spreadsheet MSF.xlsx\n";
			testMSY = 1;
		}

		//Delay difference
		//CW Dec 2015 - copied from RF May 22 2013
		// command line option for implementing delay difference model "-delaydiff"
		delaydiff=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-delaydiff",opt))>-1)
		{
			delaydiff=1;
			LOG<<"______________________________________________________\n"<<'\n';
			LOG<<"    **Implementing Delay Difference model** "<<'\n';
			LOG<<"______________________________________________________"<<'\n';
		}


	END_CALCS

  // Set up the datafile to be the name found on the first line in
  //  the original datafile (iscam.dat)
	!! ad_comm::change_datafile_name(DataFile);

	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS
	// |---------------------------------------------------------------------------------|
	// |
	// | area   f
	// | group  g
	// | sex    h
	// | year   i
	// | age    j
	// | gear   k  - number of gears with unique selectivity
  // They all have to be on seperate lines due to limitations of the tpl2cpp scanner.
	int f;
  int g;
  int h;
  int i;
  int j;
  int k;

	init_int narea;
	init_int ngroup;
	init_int nsex;
	init_int syr;
	init_int nyr;
	init_int sage;
	init_int nage;
	init_int ngear;
	vector   age(sage,nage);

	// |---------------------------------------------------------------------------------|
	// | LINKS TO MANAGE ARRAY INDEXING
	// |---------------------------------------------------------------------------------|
	// | - n_ags: total number of areas * groups * sex
	// | - n_ag:  total number of areas * groups
	// | - n_gs:  total number of groups * sex
	// | - n_area:  vector of indexes for area for each sex & group combination.
	// | - n_group: vector of indexes for stock for each sex & area combination.
	// | - n_sex:   vector of indexes for sex foe each area & group combination.
	// | - pntr_ag: matrix of indices for area and group.
	// | - pntr_gs: matrix of indices for group and sex.
	// | - pntr_ags: d3_array of indices for area group sex.
	int n_ags;
	!! n_ags = narea * ngroup * nsex;
	int n_ag;
	!! n_ag  = narea * ngroup;
	int n_gs;
	!! n_gs  = ngroup * nsex;
	ivector   n_area(1,n_ags);
	ivector  n_group(1,n_ags);
	ivector    n_sex(1,n_ags);
	imatrix  pntr_ag(1,narea,1,ngroup);
	imatrix  pntr_gs(1,ngroup,1,nsex);
	3darray pntr_ags(1,narea,1,ngroup,1,nsex);

	LOC_CALCS
		age.fill_seqadd(sage,1);
		int ig,ih,is;
		ig = 0;
		ih = 0;
		is = 0;
		for(f=1; f<=narea; f++){
			for(g=1; g<=ngroup; g++){
				ih ++;
				pntr_ag(f,g) = ih;
				for(h=1;h<=nsex;h++){
					ig ++;
					n_area(ig)  = f;
					n_group(ig) = g;
					n_sex(ig)   = h;
					pntr_ags(f,g,h) = ig;
					if(f==1){
						is ++;
						pntr_gs(g,h) = is;
					}
				}
			}
		}


		if(!mseFlag){
      LOG<<"| ----------------------- |\n";
		  LOG<<"| MODEL DIMENSION         |\n";
		  LOG<<"| ----------------------- |\n";
		  LOG<<"| narea  \t"<<narea<<'\n';
		  LOG<<"| ngroup \t"<<ngroup<<'\n';
		  LOG<<"| nsex   \t"<<nsex<<'\n';
		  LOG<<"| syr    \t"<<syr<<'\n';
		  LOG<<"| nyr    \t"<<nyr<<'\n';
		  LOG<<"| sage   \t"<<sage<<'\n';
		  LOG<<"| nage   \t"<<nage<<'\n';
		  LOG<<"| ngear  \t"<<ngear<<'\n';
		  LOG<<"| n_area \t"<<n_area<<'\n';
		  LOG<<"| n_group\t"<<n_group<<'\n';
		  LOG<<"| n_sex  \t"<<n_sex<<'\n';
		  LOG<<"| pntr_ag\t"<<pntr_ag<<'\n';
		  LOG<<"| pntr_gs\t"<<pntr_gs<<'\n';
		  LOG<<"| pntr_ags\t"<<pntr_ags(1)<<'\n';
		  LOG<<"| ----------------------- |\n\n";
		}

		/* Check for dimension errors in projection control file. */
		if(pf_cntrl(1)<syr || pf_cntrl(3)<syr || pf_cntrl(5)<syr){
			LOG<<"WARNING: start year in projection file control is"
			" less than initial model year. Setting to syr.\n";
			// exit(1);
			pf_cntrl(1) = syr;
			pf_cntrl(3) = syr;
			pf_cntrl(5) = syr;
		}
		if(pf_cntrl(2)>nyr || pf_cntrl(4)>nyr || pf_cntrl(6)>nyr){
			LOG<<"ERROR: last year in projection file control is"
			" greater than last model year.\n";
			exit(1);
		}


	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | Allocation for each gear in (ngear), use 0 for survey gears.
	// |---------------------------------------------------------------------------------|
	// | fsh_flag is used to determine which fleets should be in MSY-based ref points
	// | If dAllocation >0 then set fish flag =1 else 0
	// | nfleet is the number of non-survey gear fleet with dAllocations > 0
	// |
	int nfleet;
	init_vector dAllocation(1,ngear);
	//init_ivector catch_sex_composition(1,ngear);
	//init_ivector catch_type(1,ngear);

	ivector fsh_flag(1,ngear);
	LOC_CALCS
		dAllocation = dAllocation/sum(dAllocation);
		for(int k=1;k<=ngear;k++){
			if(dAllocation(k)>0)
				fsh_flag(k)=1;
			else
				fsh_flag(k)=0;
		}
		nfleet = sum(fsh_flag);
	END_CALCS

	ivector nFleetIndex(1,nfleet);
	LOC_CALCS
		j = 1;
		for(int k=1; k<=ngear;k++){
			if(fsh_flag(k)) nFleetIndex(j++) = k;
		}

	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | Growth and maturity parameters
	// |---------------------------------------------------------------------------------|
	// | n_ags -> number of areas * groups * sex
	init_vector d_linf(1,n_ags);
	init_vector d_vonbk(1,n_ags);
	init_vector d_to(1,n_ags);
	init_vector d_a(1,n_ags);
	init_vector d_b(1,n_ags);
	init_vector d_ah(1,n_ags);
	init_vector d_gh(1,n_ags);
	init_int n_MAT;

	int t1;
	int t2;
	LOC_CALCS

		if(n_MAT){
			t1 = sage;
			t2 = nage;
		}else{
			t1 = 0;
			t2 = 0;
		}
	END_CALCS
	init_vector d_maturityVector(t1,t2);
	matrix la(1,n_ags,sage,nage);		//length-at-age
	matrix wa(1,n_ags,sage,nage);		//weight-at-age
	matrix ma(1,n_ags,sage,nage);		//maturity-at-age
	LOC_CALCS
		if(!mseFlag){
		  LOG<<setw(8)<<setprecision(4)<<'\n';
	    LOG<<"| ----------------------- |\n";
		  LOG<<"| GROWTH PARAMETERS       |\n";
		  LOG<<"| ----------------------- |\n";
		  LOG<<"| d_linf  \t"<<d_linf<<'\n';
	  	LOG<<"| d_vonbk \t"<<d_vonbk<<'\n';
	  	LOG<<"| d_to    \t"<<d_to<<'\n';
	  	LOG<<"| d_a     \t"<<d_a<<'\n';
	  	LOG<<"| d_b     \t"<<d_b<<'\n';
	  	LOG<<"| d_ah    \t"<<d_ah<<'\n';
	  	LOG<<"| d_gh    \t"<<d_gh<<'\n';
	  	LOG<<"| ----------------------- |\n\n";
    }
	  // length & weight-at-age based on input growth pars
	  ma.initialize();
    for(ig=1;ig<=n_ags;ig++){
      la(ig) = d_linf(ig)*(1. - exp(-d_vonbk(ig)*(age-d_to(ig))));
	  	wa(ig) = d_a(ig) * pow(la(ig),d_b(ig));
	  	h = n_sex(ig);
	  	if(n_MAT==0){
	  			ma(ig) = plogis(age,d_ah(ig),d_gh(ig));
	  	}else if(n_MAT>0 && h !=2){
        ma(ig) = d_maturityVector;
	  	}
	  }


	END_CALCS

	//=======================================================================================
	//Delay difference parameters-  all fixed??
	//pergunta - discuss wit Robyn F the dimensions of these quanitites
	// RFUpdate : same as growth pars d_vonbk etc so I guess (1,n_ags)
	//=======================================================================================
	
	//age at knife-edge recruitment 
	init_ivector kage(1,n_gs); 
	
    //DD growth parameters :: RF (02-Apr-2013)
    init_vector alpha_g(1,n_gs);  //growth alpha (intercept of Ford-Walford plot; derived from wk and wk-1, H&W 1992, p 334)
	
	//growth rho (slope of Ford-Walford plot; H&W 1992, p 332)
	init_vector rho_g(1,n_gs);  

  	init_vector wk(1,n_gs);


	// |---------------------------------------------------------------------------------|
	// | Historical removal
	// |---------------------------------------------------------------------------------|
	// | - Total catch in weight (type=1), numbers (type=2), or roe (type=3).
	// | - dCatchData matrix cols: (year gear area group sex type value).
	// | - If total catch is asexual (sex=0), pool predicted catch from nsex groups.
	// | - ft_count  -> Number of estimated fishing mortality rate parameters.
	// | - d3_Ct     -> An array of observed catch in group(ig) year (row) by gear (col)
	// | - [?] - TODO: fix special case where nsex==2 and catch sex = 0 in catch array.
	init_int nCtNobs;
	init_matrix dCatchData(1,nCtNobs,1,7);
	3darray d3_Ct(1,n_ags,syr,nyr,1,ngear);

	int ft_count;

	LOC_CALCS

		ft_count = nCtNobs;

		
		if(!mseFlag){
     	  LOG<<"| ----------------------- |\n";
		    LOG<<"| HEAD(dCatchData)        |\n";
		    LOG<<"| ----------------------- |\n";
        LOG<<dCatchData.sub(1,3)<<'\n';
		    LOG<<"| ----------------------- |\n\n";
		    LOG<<"| ----------------------- |\n";
		    LOG<<"| TAIL(dCatchData)        |\n";
		    LOG<<"| ----------------------- |\n";
		    LOG<<dCatchData.sub(nCtNobs-3,nCtNobs)<<'\n';
		    LOG<<"| ----------------------- |\n";
		}

		
		
		d3_Ct.initialize();
		int k;



		for(int ii=1;ii<=nCtNobs;ii++){
			i = dCatchData(ii)(1);
			k = dCatchData(ii)(2);
			f = dCatchData(ii)(3);
			g = dCatchData(ii)(4);
			h = dCatchData(ii)(5);
			if(h==0){
        for(h=1;h<=nsex;h++){
					ig = pntr_ags(f,g,h);
					d3_Ct(ig)(i)(k) = 1./nsex*dCatchData(ii)(7);
				}
			}else{
				ig = pntr_ags(f,g,h);
				d3_Ct(ig)(i)(k) = dCatchData(ii)(7);
			}
		}

		

	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | RELATIVE ABUNDANCE INDICIES (ragged array)
	// |---------------------------------------------------------------------------------|
	// | nItNobs     = number of independent surveys
	// | n_it_nobs   = number of survey observations
	// | n_survey_type = 1: survey is proportional to vulnerable numbers
	// | n_survey_type = 2: survey is proportional to vulnerable biomass
	// | n_survey_type = 3: survey is proportional to vulnerable spawning biomass
	// | d3_survey_data: (iyr index(it) gear area group sex wt timing)
	// | it_wt       = relative weights for each relative abundance normalized to have a
	// |               mean = 1 so rho = sig^2/(sig^2+tau^2) holds true in variance pars.
	// |

	init_int nItNobs;
	ivector        nSurveyIndex(1,nItNobs);
	init_ivector      n_it_nobs(1,nItNobs);
	init_ivector  n_survey_type(1,nItNobs);
	init_3darray d3_survey_data(1,nItNobs,1,n_it_nobs,1,8);
	matrix                it_wt(1,nItNobs,1,n_it_nobs);
	matrix               it_grp(1,nItNobs,1,n_it_nobs);

	LOC_CALCS
		if(!mseFlag)
		{
		LOG<<"| ----------------------- |\n";
		LOG<<"| TAIL(d3_survey_data)    |\n";
		LOG<<"| ----------------------- |\n";
		LOG<<d3_survey_data(nItNobs).sub(n_it_nobs(nItNobs)-3,n_it_nobs(nItNobs))<<'\n';
		LOG<<"| ----------------------- |\n\n";
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

	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | AGE OR LENGTH COMPOSITION DATA (ragged object)
	// |---------------------------------------------------------------------------------|
	// | - nAgears    -> number of age-composition matrixes, one for each gear.
	// | - n_A_nobs   -> ivector for number of rows in age composition (A) matrix
	// | n_A_sage     -> imatrix for starting age in each row
	// | n_A_nage	  -> imatrix for plus group age in each row
	// | inp_nscaler  -> effective sample size for iterative re-weighting in multinomial.
	// | icol_A       -> number of columns for each row in A.
	// | A            -> array of data (year,gear,area,group,sex|Data...)
	// | d3_A_obs     -> array of catch-age data only.
	// |
	init_int nAgears;
	init_ivector n_A_nobs(1,nAgears);
	init_ivector n_A_sage(1,nAgears);
	init_ivector n_A_nage(1,nAgears);
	init_vector  inp_nscaler(1,nAgears);
	init_ivector n_ageFlag(1,nAgears);
  // The 5 in the next command is to remove the first 5 columns
  // from the age comp 'data' because they are not the actual ages,
  // but the header data.
	init_3darray d3_A(1,nAgears,1,n_A_nobs,n_A_sage-5,n_A_nage);
	3darray d3_A_obs(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	LOC_CALCS
		if( n_A_nobs(nAgears) > 0) 
		{
			if(!mseFlag)
			{
				if(n_A_nobs(nAgears) > 3)
				{
					LOG<<"| ----------------------- |\n";
					LOG<<"| TAIL(A)       |"<<'\n';
					LOG<<"| ----------------------- |\n";
					LOG<<setw(4)<<d3_A(nAgears).sub(n_A_nobs(nAgears)-2,n_A_nobs(nAgears))<<'\n';
					LOG<<"| ----------------------- |\n\n";
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
			LOG<<"| ----------------------- |\n";
			LOG<<"| NO AGE OR LENGTH DATA   |\n";
			LOG<<"| ----------------------- |\n";
		}

	END_CALCS
	// |---------------------------------------------------------------------------------|
	// | EMPIRICAL WEIGHT_AT_AGE DATA
	// |---------------------------------------------------------------------------------|
	// | Mean weight-at-age data (kg) if nWtNobs > 0
	// | sage-5 = year
	// | sage-4 = gear
	// | sage-3 = area
	// | sage-2 = stock
	// | sage-1 = sex
	// | - construct and fill weight-at-age matrix for use in the model code  (d3_wt_avg)
	// | - construct and fill weight-at-age dev matrix for length-based selex (d3_wt_dev)
	// | - construct and fill fecundity-at-age matrix for ssb calculations.   (d3_wt_mat)
	// | [ ] - TODO fix h=0 option for weight-at-age data
	// | [ ] - TODO need to accomodate ragged arrays, or NA values, or partial d3_wt_avg.
	// | [ ] - TODO Construct AgeLength TM in data section for empirical weight-at-age.
	// | nWtTab  = number of Empirical weight-at-age tables.
	// | nWtNobs = number of rows in each weight-at-age table.
	// | d3_inp_wt_avg = input weight-at-age.

	init_int nWtTab;
	init_ivector nWtNobs(1,nWtTab);
	init_3darray d3_inp_wt_avg(1,nWtTab,1,nWtNobs,sage-5,nage);
	vector tmp_nWtNobs(1,nWtTab);
	int sum_tmp_nWtNobs; 
	vector projwt(1,nWtTab);
	vector n_bf_wt_row(1,nWtTab);

	// |---------------------------------------------------------------------------------|
	// | EMPIRICAL ANNUAL MEAN WEIGHT DATA (RAGGED ARRAY)
	// |---------------------------------------------------------------------------------|
	// | FITTING TO ANNUAL MEAN WEIGHTS REQUIRES SWITCHING THIS ON AND SETTING VARIANCE  |
	// |---------------------------------------------------------------------------------|
	// | nMeanWt     - number of independent annual mean weight series
	// | nMeanWtNobs - vector :: number of annual mean weight observations in each series
  	init_int nMeanWt;
	init_ivector nMeanWtNobs(1,nMeanWt);
	init_3darray d3_mean_wt_data(1,nMeanWt,1,nMeanWtNobs,1,8);

	LOC_CALCS
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
	END_CALCS

		3darray xinp_wt_avg(1,nWtTab,1,tmp_nWtNobs,sage-5,nage);
		matrix  xxinp_wt_avg(1,sum_tmp_nWtNobs,sage-5,nage);

	LOC_CALCS

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
	END_CALCS

	matrix  dWt_bar(1,n_ags,sage,nage);
	3darray d3_wt_avg(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_dev(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_mat(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_len_age(1,n_ags,syr,nyr+1,sage,nage);

	// Trying to figure this out for Robyn forrest.
	// imatrix nrh(1,2,1,2);
	// imatrix nch(1,2,1,2);
	// !!nrh(1) = 1;
	// !!nrh(2) = 2;
	// !!nch(1) = 1;
	// !!nch(2) = 2;
	// !!LOG<<nrh;
	// !!LOG<<nch;
	// 4darray d4_alk(1,nAgears,1,n_A_nobs,sage,nrh,nage,nch);

	LOC_CALCS

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
				//LOG<<sum(first_difference(mtmp(j)(syr,nyr)));
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
				LOG<<"|-----------------------------------------------|\n";
				LOG<<"| ERROR IN INPUT DATA FILE FOR MEAN WEIGHT DATA |\n";
				LOG<<"|-----------------------------------------------|\n";
				LOG<<"| - Cannot have an observed mean weight-at-age  |\n";
				LOG<<"|   less than or equal to 0.  Please fix.       |\n";
				LOG<<"| - You are permitted to use '-99.0' for missing|\n";
				LOG<<"|   values in your weight-at-age data.          |\n";
				LOG<<"| - Aborting program!                           |\n";
				LOG<<"|-----------------------------------------------|\n";
				ad_exit(1);
			}
		}

	END_CALCS
	// |---------------------------------------------------------------------------------|
	// | END OF DATA FILE
	// |---------------------------------------------------------------------------------|
	// |
	init_int eof;
	LOC_CALCS
		
	  if(eof==999){
      	LOG<<" ______________________________ \n";
     	LOG<<"|      END OF DATA SECTION     |\n";
	    LOG<<"|______________________________|\n";
	  }else{
      	LOG<<"\n *** ERROR READING DATA *** \n\n";
      	exit(1);
	  }

	END_CALCS
	// |---------------------------------------------------------------------------------|
	// | VARIABLES FOR MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	matrix fmsy(1,ngroup,1,nfleet);	//Fishing mortality rate at Fmsy
	matrix fall(1,ngroup,1,nfleet);	//Fishing mortality based on dAllocation
	matrix  msy(1,ngroup,1,nfleet);	//Maximum sustainable yield
	vector bmsy(1,ngroup);			//Spawning biomass at MSY
 // number Umsy;					//Exploitation rate at MSY
	vector age_tau2(1,nAgears);	//MLE estimate of the variance for age comps
 // 	//catch-age for simulation model (could be declared locally 3d_array)
 // 	3darray d3C(1,ngear,syr,nyr,sage,nage);		
	
	
		
	
	// |---------------------------------------------------------------------------------|
	// | CONTROL FILE
	// |---------------------------------------------------------------------------------|
	// |
	!! ad_comm::change_datafile_name(ControlFile);
	

	// |---------------------------------------------------------------------------------|
	// | Leading Parameters
	// |---------------------------------------------------------------------------------|
	// | npar            -> number of leading parameters
	// | ipar_vector     -> integer vector based on the number of areas groups sexes
	// | -1) log_ro      - unfished sage recruitment
	// | -2) steepness   - steepness of the stock-recruitment relationship
	// | -3) log_m       - instantaneous natural mortality rate
	// | -4) log_avgrec  - average sage recruitment from syr+1 to nyr
	// | -5) log_recinit - average sage recruitment for initialization
	// | -6) rho         - proportion of total variance for observation errors
	// | -7) vartheta    - total precision (1/variance)
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	
	vector   theta_ival(1,npar);
	vector     theta_lb(1,npar);
	vector     theta_ub(1,npar);
	ivector   theta_phz(1,npar);
	ivector theta_prior(1,npar);
	ivector ipar_vector(1,npar);
	LOC_CALCS
		theta_ival  = column(theta_control,1);
		theta_lb    = column(theta_control,2);
		theta_ub    = column(theta_control,3);
		theta_phz   = ivector(column(theta_control,4));
		theta_prior = ivector(column(theta_control,5));
		ipar_vector(1,2) = ngroup;
		ipar_vector(6,7) = ngroup;
		ipar_vector(3)   = n_gs;
		ipar_vector(4,5) = n_ag;
	END_CALCS
	
	// |---------------------------------------------------------------------------------|
	// | CONTROLS PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears                  |
	// |---------------------------------------------------------------------------------|
	// |
	
	init_ivector nCompIndex(1,nAgears);
	init_ivector nCompLikelihood(1,nAgears);
	init_vector  dMinP(1,nAgears);
	init_vector  dEps(1,nAgears);
	init_ivector nPhz_age_tau2(1,nAgears);
	init_ivector nPhz_phi1(1,nAgears);
	init_ivector nPhz_phi2(1,nAgears);
	init_ivector nPhz_df(1,nAgears);
	init_int check;
	LOC_CALCS
		LOG<<"check is "<<check<<'\n';
		if(check != -12345) 
		{
			LOG<<"check is "<<check<<'\n';
			LOG<<"Check integer for EOF, should be -12345.. = "<<check<<'\n';
      		LOG<<"Error reading composition controls\n";
      		exit(1);
		}
	END_CALCS


	// |---------------------------------------------------------------------------------|
	// | CONTROLS FOR SELECTIVITY OPTIONS
	// |---------------------------------------------------------------------------------|
	// | - 12 different options for modelling selectivity which are summarized here:
	// | - isel_npar  -> ivector for # of parameters for each gear.
	// | - jsel_npar  -> ivector for the number of rows for time-varying selectivity.
	// | 
	// | SEL_TYPE  DESCRIPTION
	// |    1      age-based logistic function with 2 parameters.
	// |    2      age-based selectivity coefficients with nage-sage parameters.
	// |    3      cubic spline with age knots.
	// |    4      time-varying cubic spline with age knots.
	// |    5      time-varying bicubic spline with age and year knots.
	// |    6      logistic with fixed parameters.
	// |    7      logistic function of body weight with 2 parameters.
	// |    8      logistic 3 parameter function based on mean weight deviations.
	// |    11     length-based logistic function with 2 parametrs based on mean length.
	// |    12     length-based selectivity coefficients with cubic spline interpolation.
	// |	13 	   age-based selectivity coefficients with age_min-age_max parameters.
	// |
	// | selex_controls (1-10)
	// |  1  -> isel_type - switch for selectivity.
	// |  2  -> ahat (sel_type=1) - age-at-50% vulnerbality for logistic function or (sel_type=13) -age_min
	// |  3  -> ghat (sel_type=1) - std at 50% age of vulnerability for logistic function or (sel_type=13) -age_max
	// |  4  -> age_nodes - No. of age-nodes for bicubic spline.
	// |  5  -> yr_nodes  - No. of year-nodes for bicubic spline.
	// |  6  -> sel_phz   - phase for estimating selectivity parameters.
	// |  7  -> lambda_1  - penalty weight for 2nd difference in selectivity.
	// |  8  -> lambda_2  - penalty weight for dome-shaped selectivity.
	// |  9  -> lambda_3  - penalty weight for 2nd difference in time-varying selectivity.
	// |  10 -> Number of discrete selectivity blocks.
	// |

	init_matrix selex_controls(1,10,1,ngear);
	

	ivector    isel_npar(1,ngear);	
	ivector    jsel_npar(1,ngear);	
	ivector    isel_type(1,ngear);	
	ivector      sel_phz(1,ngear);	
	ivector n_sel_blocks(1,ngear);	
	vector      ahat_agemin(1,ngear);	
	vector      ghat_agemax(1,ngear);
	vector age_nodes(1,ngear);	
	vector  yr_nodes(1,ngear);	
	vector  lambda_1(1,ngear);	
	vector  lambda_2(1,ngear);	
	vector  lambda_3(1,ngear);	
	
	LOC_CALCS
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
	END_CALCS
	
	init_imatrix sel_blocks(1,ngear,1,n_sel_blocks);


	LOC_CALCS
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
	END_CALCS
	
 	// |---------------------------------------------------------------------------------|
	// | PRIOR FOR RELATIVE ABUNDANCE DATA
	// |---------------------------------------------------------------------------------|
	// | nits     -> number of relative abundance indices
	// | q_prior  -> type of prior to use, see legend
	// | mu_log_q -> mean q in log-space
	// | sd_log_q -> std of q prior in log-space
	// |
	// | q_prior type:
	// | 0 -> uninformative prior.
	// | 1 -> normal prior on q in log-space.
	// | 2 -> penalized random walk in q.
	init_int nits;					
	init_ivector q_prior(1,nits);
	init_vector mu_log_q(1,nits);
	init_vector sd_log_q(1,nits);

	// |---------------------------------------------------------------------------------|
	// | ## CONTROLS FOR FITTING TO MEAN WEIGHT DATA	  //START_RF_ADD
	// | ## RF added this for testing with Pacific Cod data 
	// |---------------------------------------------------------------------------------|
	// |  fitMeanWt    SWITCH :: 1= fit to annual mean weights; 0=do not fit to annual mean weights
	// |  nMeanWtCV  Number of annual mean weight series
	 // |  vector weight_sig   sd for likelihood for fitting to annual mean weight (one for each series) 
	 init_int fitMeanWt;
	 init_int nMeanWtCV;
	 init_vector weight_sig(1,nMeanWtCV);                        		       //END_RF_ADD

	// |---------------------------------------------------------------------------------|
	// | Miscellaneous controls                                                          |
	// |---------------------------------------------------------------------------------|
	// | 1 -> verbose
	// | 2 -> recruitment model (1=beverton-holt, 2=rickers)
	// | 3 -> std in catch first phase
	// | 4 -> std in catch in last phase
	// | 5 -> assumed unfished in first year (0=FALSE, 1=TRUE)
	// | 6 -> Maternal effects power parameter.
	// | 7 -> mean fishing mortality rate to regularize the solution
	// | 8 -> standard deviation of mean F penalty in first phases
	// | 9 -> standard deviation of mean F penalty in last phase.
	// | 10-> phase for estimating deviations in natural mortality.
	// | 11-> std in natural mortality deviations.
	// | 12-> number of estimated nodes for deviations in natural mortality
	// | 13-> fraction of total mortality that takes place prior to spawning
	// | 14-> number of prospective years to start estimation from syr.
	// | 15-> switch for generating selex based on IFD and cohort biomass
	init_vector d_iscamCntrl(1,16);
	int verbose;
	init_int eofc;
	LOC_CALCS
		verbose = d_iscamCntrl(1);
		if(verbose) LOG<<d_iscamCntrl;
		for(int ig=1;ig<=n_ags;ig++)
		{
			for(int i = syr; i <= nyr; i++)
			{
				d3_wt_mat(ig)(i) = pow(d3_wt_mat(ig)(i),d_iscamCntrl(6));
			}
		}

		if(eofc==999){
      LOG<<" ______________________________ \n";
	    LOG<<"|     END OF CONTROL FILE      |\n";
	    LOG<<"|______________________________|\n";
		}else{
			LOG<<"\n ***** ERROR READING CONTROL FILE ***** \n";
      exit(1);
		}
	END_CALCS

	int nf;
	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENSIONS FOR NEGATIVE LOG LIKELIHOODS   --original declaration of ilvec
	// |---------------------------------------------------------------------------------|
	// | ilvec[1,5,6,7] -> number of fishing gears (ngear)
	// | ilvec[2]       -> number of surveys       (nItNobs)
	// | ilvec[3]       -> number of age-compisition data sets (nAgears)
	// | ilvec[4]       -> container for recruitment deviations.
	//uncomment 6 lines below to revert to original and delete RF copy below
	//ivector ilvec(1,7);
	//!! ilvec    = ngear;
	//!! ilvec(1) = 1;
	//!! ilvec(2) = nItNobs;
	//!! ilvec(3) = nAgears;
	//!! ilvec(4) = ngroup;

	//RF Added extra likelihood component for annual mean weight data
	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENSIONS FOR NEGATIVE LOG LIKELIHOODS
	// |---------------------------------------------------------------------------------|
	// | ilvec[1,6,7,8] -> number of fishing gears (ngear)
	// | ilvec[2]       -> number of surveys       (nItNobs)
	// | ilvec[3]       -> number of age-compisition data sets (nAgears)
	// | ilvec[4]       -> container for recruitment deviations.
	// | ilvec[5]       -> number of annual mean weight datasets.
	ivector ilvec(1,8);
	!! ilvec    = ngear;
	!! ilvec(1) = 1;
	!! ilvec(2) = nItNobs;
	!! ilvec(3) = nAgears;
	!! ilvec(4) = ngroup;
	!! ilvec(5) = nMeanWt;

	//Different components for delay difference model
	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENSIONS FOR NEGATIVE LOG LIKELIHOODS DELDIFF
	// |---------------------------------------------------------------------------------|
	// | ilvec[1,6,7,8] -> number of fishing gears (ngear)
	// | ilvec[2]       -> number of surveys       (nItNobs)
	// | ilvec[3]       -> number of age-compisition data sets (nAgears)
	// | ilvec[4]       -> container for recruitment deviations.
	// | ilvec[5]       -> number of annual mean weight datasets.
	
		ivector ilvec_dd(1,4);
		
		!! ilvec_dd(1) = 1;
		!! ilvec_dd(2) = nItNobs;
		!! ilvec_dd(3) = ngroup;
		!! ilvec_dd(4) = nMeanWt;



	// |---------------------------------------------------------------------------------|
	// | RETROSPECTIVE ADJUSTMENT TO nyrs
	// |---------------------------------------------------------------------------------|
	// | - Don't read any more input data from here on in.
	// | - Modifying nyr to allow for retrospective analysis.
	// | - If retro_yrs > 0, then ensure that pf_cntrl arrays are not greater than nyr,
	// |   otherwise arrays for mbar will go out of bounds.
	// | - Reduce ft_count so as not to bias estimates of ft.
	// | - Establish retrospective counter for Composition data n_naa;

	ivector n_naa(1,nAgears);

	!! nyr = nyr - retro_yrs;
	LOC_CALCS
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
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | PROSPECTIVE ADJUSTMENT TO syr                                                   |
	// |---------------------------------------------------------------------------------|
	// | - start assessment at syr + # of prospective years.
	// | - adjust sel_blocks to new syr
	// | - Reduce ft_count so as not to bias estimates of ft.
	// | - Establish prospective counter for Composition data n_saa;

	ivector n_saa(1,nAgears);
	!! syr = syr + (int)d_iscamCntrl(14);

	LOC_CALCS

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
	END_CALCS

	!! LOG<<n_saa<<'\n';
	!! LOG<<n_naa<<'\n';

	// |---------------------------------------------------------------------------------|
	// | MANAGEMENT STRATEGY EVALUATION INPUTS
	// |---------------------------------------------------------------------------------|
	// |
	//LOC_CALCS
	//	ifstream ifile("Halibut2012.mse");
	//	if(ifile)
	//	{
	//		LOG<<"Vader is happy\n";;
	//		readMseInputs();
	//		exit(1);
	//	}
	//END_CALCS


	// END OF DATA_SECTION
	!! if(verbose) LOG<<"||-- END OF DATA_SECTION --||\n";

	// |--------------------------------------|
	// | Friend Class Operating Model for MSE |
	// |--------------------------------------|
	friend_class OperatingModel;

INITIALIZATION_SECTION
  theta theta_ival;
  phi1 0.01;

PARAMETER_SECTION
	// |---------------------------------------------------------------------------------|
	// | LEADING PARAMTERS
	// |---------------------------------------------------------------------------------|
	// | - Initialized in the INITIALIZATION_SECTION with theta_ival from control file.
	// | [ ] Change to init_bounded_vector_vector.
	// | theta[1] -> log_ro, or log_msy
	// | theta[2] -> steepness(h), or log_fmsy
	// | theta[3] -> log_m
	// | theta[4] -> log_avgrec
	// | theta[5] -> log_recinit
	// | theta[6] -> rho
	// | theta[7] -> vartheta
	// |
	init_bounded_vector_vector theta(1,npar,1,ipar_vector,theta_lb,theta_ub,theta_phz);

	// |---------------------------------------------------------------------------------|
	// | SELECTIVITY PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - This is a bounded matrix vector where the dimensions are defined by the 
	// | - selectivity options specified in the control file.
	// | - There are 1:ngear arrays, having jsel_npar rows and isel_npar columns.
	// | - If the user has not specified -ainp or -binp, the initial values are set
	// |   based on ahat and ghat in the control file for logistic selectivities.
	// | - Special case: if SimFlag=TRUE, then add some random noise to ahat.
	// | - NB  sel_par is in log space.
	// |
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-25.,25.,sel_phz);

	LOC_CALCS
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
					}
				}
			}
		}

	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY RATE PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - Estimate all fishing mortality rates in log-space.
	// | - If in simulation mode then initialize with F=0.1; Actual F is conditioned on
	// |   the observed catch.
	// |

	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,1);

	LOC_CALCS
		if(!SimFlag) log_ft_pars = log(0.10);
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | INITIAL AND ANNUAL RECRUITMENT
	// |---------------------------------------------------------------------------------|
	// | - Estimate single mean initial recruitment and deviations for each initial
	// |   cohort from sage+1 to nage. (Rinit + init_log_rec_devs)
	// | - Estimate mean overal recruitment and annual deviations from syr to nyr.
	// | - d_iscamCntrl(5) is a flag to initialize the model at unfished recruitment (ro),
	// |   if this is true, then do not estimate init_log_rec_devs
	// | [ ] - TODO add dev contstraint for rec_devs in calc_objective_function.

	!! int init_dev_phz = 2;
	!! if(d_iscamCntrl(5)) init_dev_phz = -1;
	init_bounded_matrix init_log_rec_devs(1,n_ag,sage+1,nage,-15.,15.,init_dev_phz);
	init_bounded_matrix log_rec_devs(1,n_ag,syr,nyr,-15.,15.,2);

	// |---------------------------------------------------------------------------------|
	// | DEVIATIONS FOR NATURAL MORTALITY BASED ON CUBIC SPLINE INTERPOLATION
	// |---------------------------------------------------------------------------------|
	// | - Estimating trends in natural mortality rates, where the user specified the 
	// |   number of knots (d_iscamCntrl(12)) and the std in M in the control file, and the phase
	// |   in which to estimate natural mortality devs (d_iscamCntrl(10)).  If the phase is neg.
	// |   then natural mortality rate deviates are not estimated and M is assumed const.
	// | - This model is implemented as a random walk, where M{t+1} = M{t} + dev.
	
	!! int m_dev_phz = -1;
	!!     m_dev_phz = d_iscamCntrl(10);
	!! int  n_m_devs = d_iscamCntrl(12);
	init_bounded_vector log_m_nodes(1,n_m_devs,-5.0,5.0,m_dev_phz);
	

	// |---------------------------------------------------------------------------------|
	// | CORRELATION COEFFICIENTS FOR AGE COMPOSITION DATA USED IN LOGISTIC NORMAL       |
	// |---------------------------------------------------------------------------------|
	// | log_age_tau2 is the variance of the composition errors.
	// | phi1 is the AR1 coefficient
	// | phi2 used in AR2 process.
	init_bounded_number_vector log_age_tau2(1,nAgears,-4.65,5.30,nPhz_age_tau2);
	init_bounded_number_vector phi1(1,nAgears,-1.0,1.0,nPhz_phi1);
	init_bounded_number_vector phi2(1,nAgears,0.0,1.0,nPhz_phi2);
	init_bounded_number_vector log_degrees_of_freedom(1,nAgears,0.70,10.0,nPhz_df);

	// |---------------------------------------------------------------------------------|
	// | AUTOCORRELATION IN RECRUITMENT DEVIATIONS                                       |
	// |---------------------------------------------------------------------------------|
	// | gamma_r: what fraction of the residual from year t-2 carries over to t-1.
	init_bounded_number gamma_r(0,1,-4);
	!!gamma_r = 0;

	// |---------------------------------------------------------------------------------|
	// | OBJECTIVE FUNCTION VALUE
	// |---------------------------------------------------------------------------------|
	// | - the value that ADMB will minimize, called objfun in iSCAM
	// |
	objective_function_value objfun;
	

    // |---------------------------------------------------------------------------------|
    // | POPULATION VARIABLES
    // |---------------------------------------------------------------------------------|
    // | - m_bar       -> Average natural mortality rate from syr to nyr.
    // |
	number m_bar;	///< Average natural mortality rate.			


	// |---------------------------------------------------------------------------------|
	// | POPULATION VECTORS
	// |---------------------------------------------------------------------------------|
    // | - ro          -> theoretical unfished age-sage recruits. 
    // | - bo          -> theoretical unfished spawning biomass (MSY-based ref point).
    // | - sbo         -> unfished spawning biomass at the time of spawning.
    // | - kappa       -> Goodyear recruitment compensation ratio K = 4h/(1-h); h=K/(4+K)
    // | - so          -> Initial slope (max R/S) of the stock-recruitment relationship.
    // | - beta        -> Density dependent term in the stock-recruitment relationship.
    // | - m           -> Instantaneous natural mortality rate by nsex
    // | - log_avgrec  -> Average sage recruitment(syr-nyr,area,group).
    // | - log_recinit -> Avg. initial recruitment for initial year cohorts(area,group).
	// | - log_m_devs  -> annual deviations in natural mortality.
	// | - q           -> conditional MLE estimates of q in It=q*Bt*exp(epsilon)
	// | - ct          -> predicted catch for each catch observation
	// | - eta         -> standardized log residual (log(obs_ct)-log(ct))/sigma_{ct}
    // | - rho         -> Proportion of total variance associated with obs error.
    // | - varphi      -> Total precision of CPUE and Recruitment deviations.
    // | - sig         -> STD of the observation errors in relative abundance data.
    // | - tau         -> STD of the process errors (recruitment deviations).
	// |
	 
	vector        ro(1,ngroup);
	vector        bo(1,ngroup);
	vector        no(1,ngroup);
	vector       sbo(1,ngroup);
	vector     kappa(1,ngroup);
	vector steepness(1,ngroup);
	vector        so(1,ngroup);
	vector      beta(1,ngroup);
	vector           m(1,n_gs);	
	vector  log_avgrec(1,n_ag);			
	vector log_recinit(1,n_ag);			
	vector          q(1,nItNobs);
	vector         ct(1,nCtNobs);
	vector        eta(1,nCtNobs);	
	vector log_m_devs(syr+1,nyr);
	vector    rho(1,ngroup);	
	vector varphi(1,ngroup);
	vector    sig(1,ngroup);	
	vector    tau(1,ngroup); 
	
	// |---------------------------------------------------------------------------------|
	// | MATRIX OBJECTS
	// |---------------------------------------------------------------------------------|
	// | - log_rt   -> age-sage recruitment for initial years and annual recruitment.
	// | - catch_df -> Catch data_frame (year,gear,area,group,sex,type,obs,pred,resid)
	// | - eta      -> log residuals between observed and predicted total catch.
	// | - nlvec    -> matrix for negative loglikelihoods.
	// | - epsilon  -> residuals for survey abundance index
	// | - it_hat   -> predicted survey index (no need to be differentiable)
	// | - qt       -> catchability coefficients (time-varying)
	// | - sbt      -> spawning stock biomass by group used in S-R relationship.
	// | - bt       -> average biomass by group used for stock projection
	// | - rt          -> predicted sage-recruits based on S-R relationship.
	// | - delta       -> residuals between estimated R and R from S-R curve (process err)
	// | - annual_mean_weight    ->  ragged matrix of estimated annual mean weights for each gear with empirical annual mean weight observations  //START_RF_ADD   END_RF_ADD    RF ADDED THIS FOR TESTING WITH PACIFIC COD DATA
	// |     
	
	matrix  log_rt(1,n_ag,syr-nage+sage,nyr);
	
	//START_RF_ADD
	//matrix   nlvec(1,7,1,ilvec);	 // original declaration
	matrix   nlvec(1,8,1,ilvec);	  //added extra component to objective function to incorporate annual mean weight data (also modified ilvec)
	matrix   nlvec_dd(1,4,1,ilvec_dd);	  //added extra component to objective function to incorporate annual mean weight data (also modified ilvec)
	
	//END_RF_ADD

	matrix epsilon(1,nItNobs,1,n_it_nobs);
	matrix  it_hat(1,nItNobs,1,n_it_nobs);
	matrix      qt(1,nItNobs,1,n_it_nobs);
	matrix     sbt(1,ngroup,syr,nyr+1);
	matrix      bt(1,ngroup,syr,nyr+1);
	matrix      rt(1,ngroup,syr+sage,nyr); 
	matrix   delta(1,ngroup,syr+sage,nyr);
	matrix   annual_mean_weight(1,nMeanWt,1,nMeanWtNobs);
	matrix   obs_annual_mean_weight(1,nMeanWt,1,nMeanWtNobs);

	// |---------------------------------------------------------------------------------|
	// | THREE DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | - ft       -> Mean fishing mortality rates for (area-sex, gear, year)
	// | F          -> Instantaneous fishing mortality rate for (group,year,age)
	// | M          -> Instantaneous natural mortality rate for (group,year,age)
	// | Z          -> Instantaneous total  mortalityr rate Z=M+F for (group,year,age)
	// | S          -> Annual survival rate exp(-Z) for (group,year,age)
	// | N          -> Numbers-at-age for (group,year+1,age)
	// | - A_hat    -> ragged matrix for predicted age-composition data.
	// | - A_nu		-> ragged matrix for age-composition residuals.
	// |   -vbt    -> //vulnerable biomass to all gears //Added by RF March 19 2015
	// |
	3darray  ft(1,n_ags,1,ngear,syr,nyr);
	3darray   F(1,n_ags,syr,nyr,sage,nage);
	3darray   M(1,n_ags,syr,nyr,sage,nage);
	3darray   Z(1,n_ags,syr,nyr,sage,nage);
	3darray   S(1,n_ags,syr,nyr,sage,nage);
	3darray   N(1,n_ags,syr,nyr+1,sage,nage);
	3darray  A_hat(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	3darray   A_nu(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	3darray   vbt(1,ngroup,1,ngear,syr,nyr+1);
			
	// //matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	// //matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	 
	// matrix Z(syr,nyr,sage,nage);
	// matrix S(syr,nyr,sage,nage);
	// matrix pit(1,nItNobs,1,n_it_nobs);			//predicted relative abundance index
	

	// 3darray Ahat(1,nAgears,1,n_A_nobs,n_A_sage-2,n_A_nage);		//predicted age proportions by gear & year
	// 3darray A_nu(1,nAgears,1,n_A_nobs,n_A_sage-2,n_A_nage);		//residuals for age proportions by gear & year
	
	// |---------------------------------------------------------------------------------|
	// | FOUR DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | log_sel    -> Selectivity for (gear, group, year, age)
	// | Chat       -> Predicted catch-age array for (gear, group, year, age)
	// | 
	
	4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
	// 4darray    Chat(1,ngear,1,n_ags,syr,nyr,sage,nage);		
	

	// |---------------------------------------------------------------------------------|
	// | SDREPORT VARIABLES AND VECTORS
	// |---------------------------------------------------------------------------------|
	// | sd_depletion -> Predicted spawning biomass depletion level bt/Bo
	// | sd_log_sbt   -> Log Spawning biomass for each group.
	// |
	sdreport_vector sd_depletion(1,ngroup);	
	sdreport_matrix sd_log_sbt(1,ngroup,syr,nyr+1);
	

	// |---------------------------------------------------------------------------------|
	// | delay difference model quantities
	// |---------------------------------------------------------------------------------|
	// | numbers -> numbers of fish in the delay difference model
	// | biomass -> biomass in the delay difference model
	// | surv -> //delaydiff only

	vector snat(1,n_gs); 	//natural survival
	vector sfished(1,n_ags); 	//natural survival
  	vector wbar(1,n_gs); 
	matrix numbers(1,n_ags,syr,nyr+1);
	matrix biomass(1,n_ags,syr,nyr+1);//RF added biomass in the delay difference model - in the ASM this is set to spawning biomass
	matrix vbcom(1,n_ags,syr,nyr); 		//RF added vulnerable biomass in gear 1 - commercial fishery
	matrix vncom(1,n_ags,syr,nyr); //RF added vulnerable numbers in gear 1 - commercial fishery
	matrix annual_mean_wt(1,n_ags,syr,nyr);  //RF addition for P cod 
	matrix   F_dd(1,n_ags,syr,nyr);
	matrix   M_dd(1,n_ags,syr,nyr);
	matrix   Z_dd(1,n_ags,syr,nyr);
	matrix  surv(1,n_ags,syr,nyr);	     //survival rate in each ags by age and year
	

PRELIMINARY_CALCS_SECTION
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

	if(verbose) LOG<<"||-- END OF PRELIMINARY_CALCS_SECTION --||\n";

RUNTIME_SECTION
    maximum_function_evaluations 2000,  2000,   2000, 25000, 25000
    convergence_criteria        0.01, 0.01, 1.e-3, 1.e-4, 1.e-5


PROCEDURE_SECTION
	
	if(!delaydiff){	
		if(d_iscamCntrl(5)==2) d_iscamCntrl(5)=0; //This control determines whether population is unfished in syr (0=false). The delay diff model also has option 2 where the population is at equilibrium with fishing mortality - not implemented in ASM.
	

		initParameters();
		calcSelectivities(isel_type);
		calcTotalMortality();
		calcNumbersAtAge();
		calcTotalCatch();
		calcComposition();
		calcSurveyObservations();
		calcStockRecruitment();
		calcAnnualMeanWeight();
	
	}	

	if(delaydiff){
		initParameters();
		calcTotalMortality_deldiff();
		calcNumbersBiomass_deldiff();
		calcFisheryObservations_deldiff();
		calcSurveyObservations_deldiff();
		calcStockRecruitment_deldiff();
		calcAnnualMeanWeight_deldiff(); //RF added this for P cod - only gets added to objective function if cntrl(15)==1
	}

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
    LOG<<"Running mceval phase\n";
	}
	if(verbose){
    LOG<<"End of main function calls\n";
  }

FUNCTION void calcSdreportVariables()
  {
	sd_depletion.initialize();
	sd_log_sbt.initialize();

	for(g=1;g<=ngroup;g++)
	{
		sd_depletion(g) = sbt(g)(nyr)/sbo(g);

		sd_log_sbt(g) = log(sbt(g));
	}
	if(verbose){
    LOG<<"**** Ok after calcSdreportVariables ****\n";
  }
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
FUNCTION void initParameters()
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
	
	if(verbose){
    LOG<<"**** Ok after initParameters ****\n";
  }
	
  }
	
FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
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
		LOG<<spline_nodes<<'\n';
		vcubic_spline_function test_ffa(ia,spline_nodes);
		LOG<<test_ffa(fa)<<'\n';
		exit(1);*/
	return(ffa(fa));
  }

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs, const dvector& la)
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
  
FUNCTION dvector cubic_spline(const dvector& spline_coffs, const dvector& la)
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

// FUNCTION dvar_matrix cubic_spline_matrix(const dvar_matrix& spline_coffs)
//   {
// 	RETURN_ARRAYS_INCREMENT();
// 	int nodes= spline_coffs.colmax()-spline_coffs.colmin()+1;
// 	int rmin = spline_coffs.rowmin();
// 	int rmax = spline_coffs.rowmax();
	
// 	dvector ia(1,nodes);
// 	dvector fa(sage,nage);
// 	ia.fill_seqadd(0,1./(nodes-1));
// 	//fa.fill_seqadd(sage,1);
// 	fa.fill_seqadd(0,1./(nage-sage));
// 	vcubic_spline_function_array fna(rmin,rmax,ia,spline_coffs);
// 	RETURN_ARRAYS_DECREMENT();
// 	return(fna(fa));
	
//   }


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
FUNCTION void calcSelectivities(const ivector& isel_type)
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

						// LOG<<"Testing selex class"<<'\n';
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
 
	if(verbose){
    LOG<<"**** Ok after calcSelectivities ****\n";
  }
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
FUNCTION calcTotalMortality
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
	 
	if(verbose){
    LOG<<"**** OK after calcTotalMortality ****\n";
  }
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
FUNCTION calcNumbersAtAge
  {
	int ig,ih, kgear;
	N.initialize();
	bt.initialize();
	vbt.initialize();     //vulnerable biomass to all gears //Added by RF March 19 2015
		
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
			//bt(g)(i) += N(ig)(i) * d3_wt_avg(ig)(i);
			bt(g)(i) = sum(elem_prod(N(ig)(i),d3_wt_avg(ig)(i)));  //RF changed this. CG said it's easier to understand than the overloaded += operator.
			//vulnerable biomass to all gears //Added by RF March 19 2015
			for(kgear=1; kgear<=ngear; kgear++) vbt(g)(kgear)(i) = sum(elem_prod(elem_prod(N(ig)(i),d3_wt_avg(ig)(i)), mfexp(log_sel(kgear)(ig)(i)))); 
		}
		N(ig)(nyr+1,sage) = 1./nsex * mfexp( log_avgrec(ih));	 //No deviation
		//bt(g)(nyr+1) += N(ig)(nyr+1) * d3_wt_avg(ig)(nyr+1);
		bt(g)(nyr+1) = sum(elem_prod(N(ig)(nyr+1),d3_wt_avg(ig)(nyr+1)));
		//vulnerable biomass to all gears //Added by RF March 19 2015
		for(kgear=1; kgear<=ngear; kgear++) vbt(g)(kgear)(nyr+1) = sum(elem_prod(elem_prod(N(ig)(nyr+1),d3_wt_avg(ig)(nyr+1)), mfexp(log_sel(kgear)(ig)(nyr))));      //use nyr selectivity
	}
	if(verbose)LOG<<"**** Ok after calcNumbersAtAge ****\n";
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
  	
FUNCTION calcComposition
  {
  	int ii,ig,kk;
    ig = 0;
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
          fa = va;
				}
				else
				{
					fa = ft(ig)(k)(i) * va;
				}
        ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);
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
            fa = va;
					}
					else
					{
						fa = ft(ig)(k)(i) * va;
					}
          ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);
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

	if(verbose){
    LOG<<"**** Ok after calcComposition ****\n";
  }

  }	

FUNCTION calcTotalCatch
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
	if(verbose){
    LOG<<"**** Ok after calcTotalCatch ****\n";
  }
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
FUNCTION calcSurveyObservations
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
	
	if(verbose){
    LOG<<"**** Ok after calcSurveyObservations ****\n";
  }
	
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
		Pseudocode:
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
FUNCTION void calcStockRecruitment()
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
  if(verbose){
    LOG<<"**** Ok after calcStockRecruitment ****\n";
  }
	
  }



  //RF added for comparison to Pacific Cod delay difference model	      //START_RF_ADD
FUNCTION calcAnnualMeanWeight
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
	if(verbose){
    LOG<<"**** Ok after calcAnnualMeanWeight ****\n";
  }
  }

  //========================================================
  //Delay difference functions
  


FUNCTION calcTotalMortality_deldiff
  {

  	/**
  	Purpose: This function calculates fishing mortality, total mortality and annual
  	         surivival rates S=exp(-Z) for each year  and area*sex*group.  Z also is updated with time-varying 
  	         natural mortality rates if specificed by user.
  	Author: Catarina Wor - Adapted from Steven Martell and Robyn Forrest work
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/
  	
	int ig,ii,i,k;
	int ft_counter = 0;
	dvariable ftmp;
	
	
	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY
	// |---------------------------------------------------------------------------------|
	
	F_dd.initialize(); 

    for(ii=1;ii<=nCtNobs;ii++)
	{
		i  = dCatchData(ii)(1);	//year
		k  = dCatchData(ii)(2);  //gear
		f  = dCatchData(ii)(3);  //area
		g  = dCatchData(ii)(4);  //group
		h  = dCatchData(ii)(5);  //sex

		if( i < syr ) continue;
		if( i > nyr ) continue;
		
		ft_counter ++;
		
		if( h )
		{
			ig = pntr_ags(f,g,h);  
			ftmp = mfexp(log_ft_pars(ft_counter));	
			ft(ig)(k,i) = ftmp;
			F_dd(ig)(i) += ftmp;
		}
		else if( !h ) // h=0 case for asexual catch	       RFUpdate: are these loops backwards?
		{
			//sum over sex
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);    
				ftmp = mfexp(log_ft_pars(ft_counter));
				ft(ig)(k,i) = ftmp;
				F_dd(ig)(i) += ftmp;
					
			}
		}
	}
    

	// |---------------------------------------------------------------------------------|
	// | NATURAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// | - uses cubic spline to interpolate time-varying natural mortality
	M_dd.initialize();
	log_m_devs.initialize();
	Z_dd.initialize();
	
	int gg,hh;

	for(int igg=1;igg<=n_ags;igg++)
	{
		gg = n_group(igg);
		hh = n_sex(igg);
		M_dd(igg)(syr) = m( pntr_gs(gg,hh) );
		Z_dd(igg)(syr) = M_dd(igg)(syr) + F_dd(igg)(syr);
		surv(igg,syr) = mfexp(-Z_dd(igg,syr));

				
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
	 	 
		

		for(int iyr=syr+1; iyr<=nyr; iyr++)
		{
			M_dd(igg)(iyr) = M_dd(igg)(iyr-1); //* mfexp(log_m_devs(iyr)); RFUpdate: The delay difference model assumes constant M
			Z_dd(igg)(iyr) = M_dd(igg)(iyr) + F_dd(igg)(iyr);
			surv(igg,iyr) = mfexp(-Z_dd(igg,iyr));

		}	
		

		// |---------------------------------------------------------------------------------|
		// | TOTAL MORTALITY
		// |---------------------------------------------------------------------------------|
		// |
		
		// TODO fix for reference point calculations
		// m_bar = mean( M_tot.sub(pf_cntrl(1),pf_cntrl(2)) );	      
	
	}
	if(verbose){
    LOG<<"**** OK after  delay diff calcTotalMortality ****\n";
  }

  /*
  LOG<<"M_dd is "<<M_dd<<'\n';
  LOG<<"F_dd is "<<F_dd<<'\n';
  LOG<<"Z_dd is "<<Z_dd<<'\n';
  LOG<<"surv is "<<surv<<'\n';
  */
  //LOG<<"**** OK after  delay diff calcTotalMortality"<<'\n';
  //exit(1);
  }
	
FUNCTION calcNumbersBiomass_deldiff
    {
    	/**
  	Purpose: This function calculates  total biomass and total numbers according to the delay differnce 
  	model equations from Hilborn and Walters. Quantities are calculated for each year and area*sex*group
  	Author: Catarina Wor - Adapted from Robyn Forrest work.
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/

  	int g,h, gs;
  	//int j;
  	numbers.initialize();
  	biomass.initialize();
  	sbt.initialize();      //RFUpdate added this -- was adding to previous iteration of sbt each time
  	
  	sfished.initialize();
  	snat.initialize();
  	wbar.initialize();
        
  	//DD initialization
        //Equilibrium mean weight - obtained from weq = Beq/Neq and solving for weq
	//i.e. weq = [surv(alpha.Neq + rho.Beq + wk.R] / [surv.Neq + R]
	// with substitutions Neq = Beq/weq and R = Neq(1 - surv)
	//From SJDM, also used by Sinclair in 2005 p cod assessment
	
	//take the average of natural mortality for both sexes
	snat = mfexp(-m);
	wbar= elem_div(elem_prod(snat,alpha_g)+elem_prod(wk,(1.-snat)),(1.-elem_prod(rho_g,snat)));	

	// need to change this by multiplying quantities by 0.5 instead of taking the mean
	no.initialize();
	bo.initialize();  //RFUpdate added this

	// average sexes for unfished quantities
	for(g=1; g<=ngroup; g++)
	{
		for(h=1; h<=nsex;h++)
		{
			gs = pntr_gs(g,h);
			//H&W 1992 p339
			no(g) +=  (ro(g)*1./nsex)/(1.-snat(gs));
			bo(g) +=  (ro(g)*1./nsex)/(1.-snat(gs)) * wbar(gs);		
		}

		so(g) 	  = kappa(g)*(ro(g)/bo(g));  

		switch(int(d_iscamCntrl(2)))
		{
			case 1:  // | Beverton Holt model
				beta(g)   = (kappa(g)-1.)/bo(g);				
			break;
			
			case 2:  // | Ricker model
				beta(g)   = log(kappa(g))/bo(g);
			break;
		}	
	}
		
	/*
	LOG<<"m snat wbar ro no bo"<<'\n';
        LOG<<m<<'\n';
	LOG<<snat<<'\n';
	LOG<<wbar<<'\n';
	LOG<<ro<<'\n';
	LOG<<no<<'\n';
        LOG<<bo<<'\n'; exit(1);  
       */
       
	//recruitment for projection year
	dvar_vector rnplus=mfexp(log_avgrec); //assume recruits nyr+1 average - same as for ASM
              
	int ih,ig;	
	for(ig=1;ig<=n_ags;ig++)
	{
		f  = n_area(ig);
		g  = n_group(ig);
		h  = n_sex(ig);
		ih = pntr_ag(f,g);
		gs = pntr_gs(g,h);
	        dvar_vector tmp_N(sage,nage);
	        tmp_N.initialize();
						
		switch(int(d_iscamCntrl(5)))
		{
			case 0: //Unfished and not at equilibrium - Initialise as for ASM
			        tmp_N(sage)  = mfexp( log_avgrec(ih)+log_rec_devs(ih)(syr));
				for(int j=sage+1;j<=nage;j++)
				{
					tmp_N(j)=mfexp(log_recinit(ih)+init_log_rec_devs(ih)(j))*mfexp(-M_dd(ig)(syr)*(j-sage));
				}
		        	tmp_N(nage)/=(1.-mfexp(-M_dd(ig)(syr)));
		        	numbers(ig,syr) = sum(tmp_N)* 1./nsex;
		        	log_rt(ih,syr) = log_avgrec(ih)+log_rec_devs(ih,syr); 
				//RFUpdate Correction: below biomass is the sum of weight at age x numbers at age not wbar
				biomass(ig,syr) = sum(elem_prod(tmp_N,d3_wt_avg(ig)(syr))); 
				annual_mean_wt(ig,syr) = biomass(ig,syr)/numbers(ig,syr); //wbar(gs);

			break;
		
			case 1: //start at equilibrium unfished  //check these two options are the same in the absence of fishing mortality
                                //RFUpdate Correction: No need for age structure here
				numbers(ig,syr)= no(g);
				numbers(ig,i)/=nsex;
				biomass(ig,syr) = bo(g);
				//annual_mean_wt(ig,syr) = wbar(gs);
				// corrected this to follow RF correction for case
				annual_mean_wt(ig,syr) = biomass(ig,syr)/numbers(ig,syr);

		        log_rt(ih,syr) = log(ro(g));
			
			break;	
 
			case 2: //start at equilibrium with fishing mortality - delay difference. CHECK THIS
	  		  	sfished(ig) = surv(ig,syr); //equilibrium survivorship at initial fishing mortality (gear 1 commercial fishery)
	   		  	annual_mean_wt(ig,syr) = (sfished(ig)*alpha_g(gs) + wk(gs)*(1.-sfished(ig)))/(1-rho_g(gs)*sfished(ig));
	   		  								
	   		  	biomass(ig,syr) = -(annual_mean_wt(ig,syr)*(wk(gs)*so(g)-1)+sfished(ig)*(alpha_g(gs)+
	   		  						rho_g(gs)* annual_mean_wt(ig,syr)))/
	   		  						(beta(g)*(sfished(ig)*alpha_g(gs)+sfished(ig)*rho_g(gs)* annual_mean_wt(ig,syr)- 
	   		  							annual_mean_wt(ig,syr)));
	   		  	numbers(ig,syr) = biomass(ig,syr)/annual_mean_wt(ig,syr);
	   		  	numbers(ig,i)/=nsex;
	   		  	
	   		  	//pergunta: where does the biomass eq comes from?
	   		  	// log rt originally missing from this option
	   		  	// chose log_avgrec as placeholder-- dangerous if fishing in first year and before was very high.
	   		  	log_rt(ih,syr) = log_avgrec(ih);
	   	 	   	  		
			break;

		}

		sbt(g,syr) += biomass(ig,syr);

		for(i=syr+1;i<=nyr;i++){
			
			log_rt(ih,i)=log_avgrec(ih)+log_rec_devs(ih,i); 

			//Update biomass and numbers	
			//RFUpdate Correction: don't divide both numbers and biomass by nsex
			biomass(ig,i) =surv(ig,i-1)*(rho_g(gs)*biomass(ig,i-1)+alpha_g(gs)*numbers(ig,i-1))+
								wk(gs)*mfexp(log_rt(ih,i)); // eq. 9.2.5 in HW
			numbers(ig,i)=surv(ig,i-1)*numbers(ig,i-1)+mfexp(log_rt(ih,i)); 
			numbers(ig,i)/=nsex;
			annual_mean_wt(ig,i)=biomass(ig,i)/numbers(ig,i);		//calculate predicted weight in dynamics - possible option to fit to it
				sbt(g,i) += biomass(ig,i);
		}	
	  	  //RF doesn't like this projection step - prefers to stick to projection in projection model - this one calculates recruitment inconsistently with projection model
	  	
	  	biomass(ig,nyr+1)=(surv(ig,nyr)*(rho_g(gs)*biomass(ig,nyr)+alpha_g(gs)*numbers(ig,nyr))+wk(gs)*rnplus(ih)/nsex); 
		numbers(ig,nyr+1)=surv(ig,nyr)*numbers(ig,nyr)+rnplus(ih)/nsex;
	  	
	  	sbt(g,nyr+1) += biomass(ig,nyr+1); //set spawning biomass to biomass
	}
	
	//RFUpdate - added sbo for delay diff model
	 sbo(g)=bo(g);
	
	if(verbose){
    LOG<<"**** Ok after calcNumbersBiomass_deldiff ****\n";
  	}

  	    /*
  	LOG<<"surv is "<<surv<<'\n';
  	LOG<<"biomass is "<<biomass<<'\n';
	LOG<<"numbers is "<<numbers<<'\n';
	LOG<<"sbt is "<<sbt<<'\n';
	LOG<<"log_rt is "<<log_rt<<'\n';
	LOG<<"log_rec_devs is "<<log_rec_devs<<'\n';
	LOG<<"**** Ok after calcNumbersBiomass_deldiff ****"<<'\n';
	//exit(1);
	*/
 	
  }


		
FUNCTION calcFisheryObservations_deldiff
	{

		/**
  	Purpose: This function calculates commertial catches for each year, gear and area*group*sex.
  	Author: Catarina Wor - Adapted from Steven Martell and Robyn Forrest work
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/

		int i,k,f,g,h,l,ig,ii;

		ct.initialize();
		eta.initialize();
		double d_ct;

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
					ct(ii) = biomass(ig,i)*(1-mfexp(-M_dd(ig,i)-ft(ig)(k)(i)))*(ft(ig)(k)(i)/(M_dd(ig,i) + ft(ig)(k)(i)));
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						ct(ii) += biomass(ig,i)*(1-mfexp(-M_dd(ig,i)-ft(ig)(k)(i)))*(ft(ig)(k)(i)/(M_dd(ig,i) + ft(ig)(k)(i)));						
					}
				}
			break;

			case 2:	//catch in numbers
				if( h )
				{
					ig     = pntr_ags(f,g,h);
					ct(ii) = numbers(ig,i)*(1-mfexp(-M_dd(ig,i)-ft(ig)(k)(i)))*(ft(ig)(k)(i)/(M_dd(ig,i) + ft(ig)(k)(i)));
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						ct(ii) += numbers(ig,i)*(1-mfexp(-M_dd(ig,i)-ft(ig)(k)(i)))*(ft(ig)(k)(i)/(M_dd(ig,i) + ft(ig)(k)(i)));						
					}
				}
					
			break;
			
			case 3:	//roe - NOT IMPLEMENTED FOR DELAY DIFFERENCE MODEL
						LOG<<"WARNING: CATCH TYPE 3 (ROE FISHERY) NOT IMPLEMENTED FOR DELAY DIFFERENCE MODEL:"<<'\n';
						LOG<<"USE THE AGE-STRUCTURED MODEL - TERMINATING PROGRAM"<<'\n'; exit(1);
			break;
			
		}

		// | catch residual
		eta(ii) = log(d_ct+TINY) - log(ct(ii)+TINY);
		//eta(ii) = log(d_ct) - log(ct(ii));
	}

	if(verbose){
    	LOG<<"**** Ok after calcFisheryObservations_deldiff ****\n";
  	}
        
        /*
        LOG<<"eta is "<<eta<<'\n';
  	LOG<<"ft is "<<ft<<'\n';
  	LOG<<"ct is "<<ct<<'\n';
  	LOG<<"TINY is "<<TINY<<'\n';
	LOG<<"**** Ok after calcFisheryObservations_deldiff ****"<<'\n';
	exit(1);
	*/
			
	}

FUNCTION calcSurveyObservations_deldiff
	{

		/**
  	Purpose: This function calculates predicted survey observations for each year and area*sex*group.  Z
  	Author: Catarina Wor - Adapted from Steven Martell and Robyn Forrest work
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/

		int ii,kk,ig,nz;
		double di;
		dvariable ftmp;


		epsilon.initialize();
		it_hat.initialize();

		for(kk=1;kk<=nItNobs;kk++)
		{
			dvar_vector V(1,n_it_nobs(kk));
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
				//Na.initialize();
				for(h=1;h<=nsex;h++)
				{
					ig  = pntr_ags(f,g,h);

					dvariable z = ft(ig)(k)(i)+M_dd(ig,i);
					dvariable Np = numbers(ig,i) * exp( -z * di);
					dvariable Bp = biomass(ig,i) * exp( -z * di);


					switch(n_survey_type(kk))
					{
						case 1:
							V(ii) +=Np;
						break; 
						case 2:
							V(ii) += Bp;
						break;
						case 3:
							V(ii) += Bp;
						break;
					}
				}
		
			} // end of ii loop	

			dvector     it 	= trans(d3_survey_data(kk))(2)(iz,nz);
			dvector     wt 	= trans(d3_survey_data(kk))(7)(iz,nz);
		         wt 	= wt/sum(wt);
			
			dvar_vector zt 	= log(it) - log(V(iz,nz));
			dvariable 	zbar = sum(elem_prod(zt,wt));
			//dvariable 	zbar = mean(zt); RFUpdate: this old weighting may have been incorrect but check above line with CW
			q(kk) = mfexp(zbar);
			

		// | survey residuals
		epsilon(kk).sub(iz,nz) = zt - zbar;
		it_hat(kk).sub(iz,nz) = q(kk) * V(iz,nz);

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
			it_hat(kk).sub(iz,nz) = elem_prod(qt(kk)(iz,nz),V(iz,nz));
		}
	}
  
        /*
  	LOG<<"q is  "<<q<<'\n';
  	LOG<<"epsilon is  "<<epsilon<<'\n';
  	LOG<<"it_hat is  "<<it_hat<<'\n';
  	exit(1);
        */
        if(verbose){
    		LOG<<"**** Ok after calcSurveyObservations_deldiff ****\n";
  	}

     }
		
FUNCTION calcStockRecruitment_deldiff
	{
		/**
  	Purpose: This function calculates Recruitment based on the user choice of recruitment function
  	Recruitment is calculated for each year and for each group.
  	Author: Catarina Wor - Adapted from Steven Martell and Robyn Forrest work
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/
		//for the delay difference model so and beta are calculated in the function calcNumbersBiomass_deldiff
	
		int i,ig,f,g,ih; // ,h,gs; // h and gs aren't used yet

		dvar_vector tau(1,ngroup);
		dvar_matrix tmp_rt(1,ngroup,syr,nyr);	  //need to calc recruits in all years, then offset by kage yrs
		//dvar_matrix tmp_st(1,ngroup,syr,nyr);
		int iicount=0; 

		rt.initialize();
		delta.initialize();

		for(ig=1;ig<=n_ags;ig++)
		{
			f  = n_area(ig);
			g  = n_group(ig);
			h  = n_sex(ig);
			ih = pntr_ag(f,g);
			// gs = pntr_gs(g,h); // This isn't used yet
	
			tau(g) = sqrt(1.-rho(g))*varphi(g);

			//LOG<<"tau is "<<tau<<'\n';
			
			for(i=syr; i<=nyr; i++)
			{
				iicount++;
				
				switch(int(d_iscamCntrl(2)))
				{
					case 1:  // | Beverton Holt model
				
						if(iicount <=kage(g)){ 
							tmp_rt(g)(i) =  so(g)*sbt(g)(syr)/(1.+beta(g)*sbt(g)(syr));  
						}else{
							tmp_rt(g)(i) = so(g)*sbt(g)(i-kage(g))/(1.+beta(g)*sbt(g)(i-kage(g)));
						}

				break;

						case 2:  // | Ricker model
					
							if(iicount <=kage(g)){ 
								tmp_rt(g)(i) =so(g)*sbt(g)(syr)*exp(-beta(g)*sbt(g)(syr));
							}else{
								tmp_rt(g)(i) =so(g)*sbt(g)(i-kage(g))*exp(-beta(g)*sbt(g)(i-kage(g)));
							}

					break;
				}
			
			}
				
			rt(g)    += mfexp(log_rt(ih)(syr+sage,nyr));	 //
			
			  /*
			LOG<<"so"<<so<<'\n';
			LOG<<"sbt "<<sbt<<'\n';
			LOG<<"beta "<<beta<<'\n';
			LOG<<"rt "<<rt<<'\n';
			LOG<<"tmp_rt is "<<tmp_rt<<'\n';
			 */
			delta(g) = log(rt(g))-log(tmp_rt(g)(syr+sage,nyr))+0.5*tau*tau;
	
		}	
		//LOG<<"delta is "<<delta<<'\n';
		//exit(1); 
		if(verbose){
    	LOG<<"**** Ok after calc_stock_recruitment_deldiff ****\n";
  		}

		//LOG<<"rt is "<<rt<<'\n';
		//LOG<<"tmp_rt is "<<tmp_rt<<'\n';
  	
		//LOG<<"**** Ok after calc_stock_recruitment_deldiff ****"<<'\n';
		//exit(1);
	}

		
		
FUNCTION calcAnnualMeanWeight_deldiff
	{
			/**
  	Purpose: This function calculates the mean weight of the catch for each year, gear by dividing the total
  	biomass by the total numbers .
  	Author: Catarina Wor - Adapted from Robyn Forrest work  RFUpdate -- Steve's versions don't have this
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	TODO list:
  	
  	*/


		int ii,kk,ig,nz;
		double di;

		dvariable wN;
		dvariable wB;
		dvariable ws;
			

		for(kk=1;kk<=nMeanWt;kk++)   //loop through series with empirical annual mean weight data
		{

			dvar_vector Vn(1,nMeanWtNobs(kk));	      // | Vulnerable number-at-age to gear
			dvar_vector Vb(1,nMeanWtNobs(kk));	      // | Vulnerable biomass-at-age to gear
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
					ws  = mfexp( -Z_dd(ig)(i)*di );   //accounts for survey timing
					wN  = numbers(ig)(i)*ws;
					wB  = biomass(ig)(i)*ws;
					Vn(ii) += wN;  //adds sexes
					Vb(ii) += wB;//TODO: need to replace d3wtavg for something else
				}
		 
		 		annual_mean_weight(kk)(ii) = Vb(ii)/Vn(ii);
				obs_annual_mean_weight(kk)(ii)	= d3_mean_wt_data(kk)(ii)(2);	  //fill a matrix with observed annual mean weights - makes objective function calcs easier
			}	// end of ii loop
		} // end of kk loop

		if(verbose){
    	LOG<<"**** Ok after calcAnnualMeanWeight_deldiff ****\n";
  		}

  		//LOG<<"annual_mean_weight is "<<annual_mean_weight<<'\n';
		//LOG<<"obs_annual_mean_weight is "<<obs_annual_mean_weight<<'\n';
  	//	LOG<<"**** Ok after calcAnnualMeanWeight_deldiff ****"<<'\n';
		//exit(1);
	}


	
FUNCTION calcObjectiveFunction
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



// 	int i,j,k;
// 	double o=1.e-10;
	
	nlvec.initialize();
	nlvec_dd.initialize();
	
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
		if(!delaydiff){
			nlvec(1) = dnorm(eta,0.0,sig_c);
		}else{
			nlvec_dd(1) = dnorm(eta,0.0,sig_c);
		}

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
		
		if(!delaydiff){
			nlvec(2,k)=dnorm(epsilon(k),sig_it); 
		}else{
			nlvec_dd(2,k)=dnorm(epsilon(k),sig_it);  
		}
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
	if(!delaydiff){
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
            //LOG<<'\n';
            //LOG<<"        log_age_tau2: "<<log_age_tau2<<'\n';
            //LOG<<"                   k: "<<k<<'\n';
            //LOG<<"     log_age_tau2(k): "<<log_age_tau2(k)<<'\n';
            //LOG<<"exp(log_age_tau2(k)): "<<exp(log_age_tau2(k))<<'\n';
            //LOG<<"             phi2(k): "<<phi2(k)<<'\n';
            //LOG<<" cLN_Age(expk,phi2k): "<<cLN_Age(exp(log_age_tau2(k)))<<'\n'<<'\n';
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
	}//end if(!delaydiff)
	
	
	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT LIKELIHOOD COMPONENT
	// |---------------------------------------------------------------------------------|
	// | - tau is the process error standard deviation.
	if( active(theta(1)) || active(theta(2)) )
	{
		for(g=1;g<=ngroup;g++)
		{
			if(!delaydiff){
				nlvec(4,g) = dnorm(delta(g),tau(g));
			}else{
				nlvec_dd(3,g) = dnorm(delta(g),tau(g));		
			}
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
	
	if(!delaydiff){
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
	}//end if(!delaydiff)
	//START_RF_ADD
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR ANNUAL MEAN WEIGHT DATA
	// |---------------------------------------------------------------------------------|
	// | - sig_it     -> vector of standard deviations based on relative wt for survey.
	// |  init_3darray d3_mean_wt_data(1,nMeanWt,1,nMeanWtNobs,1,7)
  if(fitMeanWt){
	  for(k=1;k<=nMeanWt;k++){
		  dvar_vector epsilon_wt = log(annual_mean_weight(k)) - log(obs_annual_mean_weight(k));
		  if(!delaydiff){
		  	nlvec(8,k) = dnorm(epsilon_wt,weight_sig(k)); //fit to annual mean weight if fitMeanWt is switched on in the control file
		  }else{
		  	nlvec_dd(4,k) = dnorm(epsilon_wt,weight_sig(k)); //fit to annual mean weight if fitMeanWt is switched on in the control file
		  }
	  }
  }
  	//LOG<<"nlvec "<<nlvec<<'\n';
  	//LOG<<"end of nlvec "<<'\n';

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
// 	//301 787 0241  Richard Methot cell phone.
// 	//ibis charles du gaulle at
// 	//01 49 19 19 20
	

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
		LOG<<nlvec<<'\n';
		LOG<<lvec<<'\n';
		LOG<<priors<<'\n';
		LOG<<pvec<<'\n';
		LOG<<qvec<<'\n';
	}
	// LOG<<nlvec;
	

	switch(delaydiff){
				case 0:		
					objfun  = sum(nlvec);
					objfun += sum(lvec);
					objfun += sum(priors);
					objfun += sum(pvec);
					objfun += sum(qvec);

				break;
					
				case 1:
					objfun  = sum(nlvec_dd);
					objfun += sum(priors);
					objfun += sum(pvec);
					objfun += sum(qvec);

					/*
					LOG<<"nlvec_dd  "<<'\n'<<nlvec_dd<<'\n';
					LOG<<"priors  "<<priors<<'\n';
					LOG<<"pvec  "<<pvec<<'\n';
					LOG<<"qvec  "<<qvec<<'\n';
					LOG<<"objfun  "<<objfun<<'\n';
					 */
				break;
					
				}
		
	nf++;

	if(verbose){
    LOG<<"**** Ok after calcObjectiveFunction ****\n";
  }
	//LOG<<"nlvec3 is "<<nlvec(3)<<'\n';
	 //LOG<<"**** Ok after calcObjectiveFunction ****"<<'\n';
  	
	// if(last_phase()){
	// 	ad_exit(1); //para!
	// }
  	
  }

// FUNCTION void equilibrium(const double& fe, const dvector& ak, const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dmatrix& va,double& re,double& ye,double& be,double& ve,double& dye_df,double& d2ye_df2)//,double& phiq,double& dphiq_df, double& dre_df)
//   {
// 	/*
// 	Equilibrium age-structured model used to determin Fmsy and MSY based reference points.
// 	Author: Steven Martell
	
// 	Comments: 
// 	This code uses a numerical approach to determine a vector of fe_multipliers
// 	to ensure that the dAllocation is met for each gear type.
	
// 	args:
// 	fe	-steady state fishing mortality
// 	ak	-dAllocation of total ye to gear k.
// 	ro	-unfished sage recruits
// 	kap	-recruitment compensation ration
// 	m	-instantaneous natural mortality rate
// 	age	-vector of ages
// 	wa	-mean weight at age
// 	fa	-mean fecundity at age
// 	va	-mean vulnerablity at age for fe gear.

	
// 	Modified args:
// 	re	-steady state recruitment
// 	ye	-steady state yield
// 	be	-steady state spawning biomass
// 	phiq		-per recruit yield
// 	dre_df		-partial of recruitment wrt fe
// 	dphiq_df	-partial of per recruit yield wrt fe
	
// 	LUCIE'S RULE: the derivative of a sum is the sum of its derivatives.
// 	Lucie says: there is some nasty calculus in here daddy, are you sure
// 	you've got it right?
	
// 	I've got it pretty close. 
	
// 	DEPRECATE THIS FUNCTION.  NOW DONE IN THE MSY CLASS
	
// 	*/
// 	int i,j,k;
// 	int nage    = max(age);
// 	int sage    = min(age);
// 	double  dre_df;
// 	double  phif;
// 	dvector lx(sage,nage);
// 	dvector lz(sage,nage);
// 	dvector lambda(1,ngear);        //F-multiplier
// 	dvector phix(1,ngear);
// 	dvector phiq(1,ngear);
// 	dvector dphiq_df(1,ngear);
// 	dvector dyek_df(1,ngear);
// 	dvector d2yek_df2(1,ngear);
// 	dvector yek(1,ngear);
// 	dmatrix qa(1,ngear,sage,nage);
// 	dmatrix xa(1,ngear,sage,nage);  //vulnerable numbers per recruit
	
// 	lx          = pow(exp(-m),age-double(sage));
// 	lx(nage)   /=(1.-exp(-m));
// 	double phie = lx*fa;		// eggs per recruit
// 	double so   = kap/phie;
// 	double beta = (kap-1.)/(ro*phie);
// 	lambda      = ak/mean(ak);	// multiplier for fe for each gear
	
	
// 	/* Must iteratively solve for f-multilier */
// 	for(int iter=1;iter<=30;iter++)
// 	{
// 		/* Survivorship under fished conditions */
// 		lz(sage)    = 1.0;
// 		lambda     /= mean(lambda);
// 		dvector fk  = fe*lambda;
// 		dvector ra  = lambda*va;
// 		dvector za  = m + fe*ra;
// 		dvector sa  = mfexp(-za);
// 		dvector oa  = 1.0 - sa;
		
		
// 		for(k=1;k<=ngear;k++)
// 		{
// 			qa(k) = elem_prod(elem_div(lambda(k)*va(k),za),oa);
// 			xa(k) = elem_prod(elem_div(va(k),za),oa);
// 		}
		
// 		double dlz_df = 0, dphif_df = 0;
// 		dphiq_df.initialize();
// 		dre_df   = 0;
// 		for(j=sage;j<=nage;j++)
// 		{
// 			if(j>sage) lz(j)  = lz(j-1) * sa(j-1);
// 			if(j>sage) dlz_df = dlz_df  * sa(j-1) - lz(j-1)*ra(j-1)*sa(j-1);
			
// 			if(j==nage)
// 			{
// 				lz(j)  = lz(j) / oa(j);
				
// 				double t4 = (-ra(j-1)+ra(j))*sa(j)+ra(j-1);
// 				dlz_df = dlz_df/oa(j) - (lz(j-1)*sa(j-1)*t4) / square(oa(j));
// 			}
			
// 			dphif_df   = dphif_df+fa(j)*dlz_df;
// 			for(k=1;k<=ngear;k++)
// 			{
// 				double t1   = lambda(k) * wa(j) *va(k,j) * ra(j) * lz(j);
// 				double t3   = -1. + (1.+za(j)) * sa(j);
// 				double t9   = square(za(j));
// 				dphiq_df(k)+= wa(j)*qa(k,j)*dlz_df + t1 * t3 / t9; 
// 			}
// 		} 
		
// 		phif   = elem_prod(lz,exp(-za*d_iscamCntrl(13)))*fa;
// 		re     = ro*(kap-phie/phif)/(kap-1.);
// 		dre_df = ro/(kap-1.0)*phie/square(phif)*dphif_df;
		
// 		/* Equilibrium yield */
// 		for(k=1;k<=ngear;k++)
// 		{
// 			phix(k)      = sum(elem_prod(elem_prod(lz,wa),xa(k)));
// 			phiq(k)      = sum(elem_prod(elem_prod(lz,wa),qa(k)));
// 			yek(k)       = fe*re*phiq(k);
// 			dyek_df(k)   = re*phiq(k) + fe*phiq(k)*dre_df + fe*re*dphiq_df(k);
// 			d2yek_df2(k) = phiq(k)*dre_df + re*dphiq_df(k);
// 		}
		
// 		/* Iterative soln for lambda */
// 		dvector pk = yek/sum(yek);
// 		dvector t1 = elem_div(ak,pk+1.e-30);
// 		lambda     = elem_prod(lambda,t1);
// 		if(abs(sum(ak-pk))<1.e-6) break;
// 	} // end of iter
// 	ve       = re*sum(elem_prod(ak,phix));
// 	be       = re*phif;
// 	ye       = sum(yek);
// 	dye_df   = sum(dyek_df);
// 	d2ye_df2 = sum(d2yek_df2);

// 	// LOG<<"EQUILIBRIUM CODE "<<setprecision(4)<<setw(2)<<fe<<setw(3)<<" "
// 	// <<ye<<setw(5)<<" "<<dye_df<<"  "<<dyek_df(1,3)<<'\n';
//   }



	
// FUNCTION void equilibrium(const double& fe,const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dvector& va,double& re,double& ye,double& be,double& phiq,double& dphiq_df, double& dre_df)
//   {
// 	/*
// 	This is the equilibrium age-structured model that is 
// 	conditioned on fe (the steady state fishing mortality rate).
	
// 	In the case of multiple fisheries, fe is to be considered as the
// 	total fishing mortality rate and each fleet is given a specified
// 	dAllocation based on its selectivity curve.  The dAllocation to 
// 	each fleet must be specified a priori.
	
// 	args:
// 	fe	-steady state fishing mortality
// 	ro	-unfished sage recruits
// 	kap	-recruitment compensation ration
// 	m	-instantaneous natural mortality rate
// 	age	-vector of ages
// 	wa	-mean weight at age
// 	fa	-mean fecundity at age
// 	va	-mean vulnerablity at age for fe gear.
// 	ak	-dAllocation of total ye to gear k.
	
// 	Modified args:
// 	re	-steady state recruitment
// 	ye	-steady state yield
// 	be	-steady state spawning biomass
// 	phiq		-per recruit yield
// 	dre_df		-partial of recruitment wrt fe
// 	dphiq_df	-partial of per recruit yield wrt fe
	
// 	FIXME add Ricker model to reference points calculations.
// 	FIXME partial derivatives for dphif_df need to be fixed when d_iscamCntrl(13)>0.
// 	*/
// 	int i;
	
// 	int nage=max(age);
// 	int sage=min(age);
// 	dvector lx=pow(exp(-m),age-double(sage));
// 	lx(nage)/=(1.-exp(-m));
// 	dvector lz=lx;
// 	dvector za=m+fe*va;
// 	dvector sa=1.-exp(-za);
// 	dvector qa=elem_prod(elem_div(va,za),sa);
	
// 	double phie = lx*fa;		//eggs per recruit
// 	double so = kap/phie;
// 	double beta = (kap-1.)/(ro*phie);
	
	
// 	double dlz_df = 0, dphif_df = 0;
// 	dphiq_df=0; dre_df=0;
// 	for(i=sage; i<=nage; i++)
// 	{
// 		if(i>sage) lz[i]=lz[i-1]*exp(-za[i-1]);
// 		if(i>sage) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1]);
// 		if(i==nage){ //6/11/2007 added plus group.
// 					lz[i]/=(1.-mfexp(-za[i]));
					
// 					dlz_df=dlz_df/(1.-mfexp(-za[i]))
// 							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
// 					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
// 				}
// 		dphif_df=dphif_df+fa(i)*dlz_df;
// 		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
// 	}
// 	//CHANGED need to account for fraction of mortality that occurs
// 	//before the spawning season in the recruitment calculation.
// 	//LOG<<"lz\t"<<elem_prod(lz,exp(-za*d_iscamCntrl(13)))<<'\n';
// 	//exit(1);
// 	//double phif = lz*fa;
// 	double phif = elem_prod(lz,exp(-za*d_iscamCntrl(13)))*fa;
// 	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
// 	re=ro*(kap-phie/phif)/(kap-1.);
// 	//LOG<<fe<<" spr ="<<phif/phie<<'\n';
// 	if(re<=0) re=0;
// 	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
// 	ye=fe*re*phiq;
// 	be=re*phif;	//spawning biomass
	
// 	//LOG<<"Equilibrium\n"<<ro<<'\n'<<re<<'\n'<<ye<<'\n';
	
//   }
	
FUNCTION void calcReferencePoints()
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
  	if(!delaydiff){

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
		

		//LOG<<"Initial Fe "<<dftry<<'\n';
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
	for( ig = 1; ig <= n_ags; ig++ ) {
    fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
		M_bar(ig)  = colsum(value(M(ig).sub(pf_cntrl(3),pf_cntrl(4))));
		M_bar(ig) /= pf_cntrl(4)-pf_cntrl(3)+1;	
	}
	for( g = 1; g <= ngroup; g++ ) {
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
		//LOG<<"group \t"<<g<<'\n';
		//exit(1);

		Msy c_msy(d_ro,d_h,M_bar,d_rho,dWt_bar,fa_bar,&d_V);
		fmsy(g) = 0.1;
		c_msy.get_fmsy(fmsy(g));
		bo = c_msy.getBo();
		bmsy(g) = c_msy.getBmsy();
		msy(g) = c_msy.getMsy();
		//LOG<<"Old Msy class\n;
    //LOG<<c_msy<<'\n';
	}
	}
	/*RF added a test of ref point calcs - runs out the model for 100 y and calculates fmsy and bmsy conditional on model parameters and data
	  Just run once in last MPD phase. Turn off after testing.*/
	//if(!mceval_phase()) run_FRP();	  //RF ran this March 18 2015 for Arrowtooth Flounder and got perfect agreement with iscam's code above
	if(delaydiff){


		LOG<<"MSY quantitied not defined for Delay difference model"<<'\n';
		if(!mceval_phase()) run_FRPdd();	  //RF ran this March 18 2015 for Arrowtooth Flounder and got perfect agreement with iscam's code above
		
	}

	if(verbose){
    	LOG<<"**** Ok after calcReferencePoints ****\n";
  	}
  }

  /**
   * This is a simple test routine for comparing the MSY class output to the 
   * MSF.xlsx spreadsheet that was used to develop the multiple fleet msy 
   * method.  Its a permanent feature of the iscam code for testing.
   */	
FUNCTION void testMSYxls()
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
	LOG<<"Initial Fe "<<dftry<<'\n';
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSY(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dfmsy = c_MSY.getFmsy(dftry);
	LOG<<"Fmsy = "<<dfmsy<<'\n';

	dvector ak(1,2);
	ak = 0.3;
	ak(2) = 1-ak(1);
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSYk(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dkmsy = c_MSYk.getFmsy(dftry,ak);
	LOG<<"Fmsy_k ="<<dkmsy<<'\n';

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
FUNCTION void simulationModel(const long& seed)
  {
	LOG<<global_parfile<<'\n';
	bool pinfile = 0;
	LOG<<"___________________________________________________\n";
	LOG<<"  **Implementing Simulation--Estimation trial**    \n";
	LOG<<"___________________________________________________\n";
	//if(norm(log_rec_devs)!=0)
	if(global_parfile)
	{
		LOG<<"\tUsing pin file for simulation\n";;
		pinfile = 1;
	}
	LOG<<"\tRandom Seed No.:\t"<< rseed<<'\n';
	LOG<<"\tNumber of retrospective years: "<<retro_yrs<<'\n';;
	LOG<<"___________________________________________________\n\n";

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

			if(!pinfile){
        LOG<<"Add stock recruitment Model\n";
      }
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
	// LOG<<d3_Ct(1)<<'\n';


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
		LOG<<d3_inp_wt_avg(1)(1)(sage-5,nage)<<'\n';

	// |---------------------------------------------------------------------------------|
	// | 11) WRITE SIMULATED DATA TO FILE
	// |---------------------------------------------------------------------------------|
	// |
	writeSimulatedDataFile();

	
// 	calcReferencePoints();
// 	//LOG<<"	OK after reference points\n"<<fmsy<<'\n';
// 	//exit(1);
// 	//	REPORT(fmsy);
// 	//	REPORT(msy);
// 	//	REPORT(bmsy);
	
	
// 	LOG<<"___________________________________________________"<<'\n';
// 	ofstream ofs("iscam.sim");
// 	ofs<<"fmsy\n"<<fmsy<<'\n';
// 	ofs<<"msy\n"<<msy<<'\n';
// 	ofs<<"bmsy\n"<<bmsy<<'\n';
// 	ofs<<"bo\n"<<bo<<'\n';
// 	ofs<<"va\n"<<va<<'\n';
// 	ofs<<"sbt\n"<<sbt<<'\n';//<<rowsum(elem_prod(N,fec))<<'\n';
// 	ofs<<"log_rec_devs\n"<<log_rec_devs<<'\n';
// 	ofs<<"rt\n"<<rt<<'\n';
// 	ofs<<"ct\n"<<obs_ct<<'\n';
// 	ofs<<"ft\n"<<trans(ft)<<'\n';
// 	//ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<'\n';
// 	ofs<<"iyr\n"<<iyr<<'\n';
// 	ofs<<"it\n"<<it<<'\n';
// 	ofs<<"N\n"<<N<<'\n';
// 	ofs<<"A\n"<<A<<'\n';
// 	ofs<<"dlog_sel\n"<<dlog_sel<<'\n';
// 	LOG<<"  -- Simuation results written to iscam.sim --\n";
// 	LOG<<"___________________________________________________"<<'\n';
	
// 	//LOG<<N<<'\n';
// 	//exit(1);
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
FUNCTION writeSimulatedDataFile
  {
  	adstring sim_datafile_name = "Simulated_Data_"+str(rseed)+".dat";
  	ofstream dfs(sim_datafile_name);
  	dfs<<"#Model dimensions"<<'\n';
  	dfs<< narea 		<<'\n';
  	dfs<< ngroup		<<'\n';
  	dfs<< nsex			<<'\n';
  	dfs<< syr   		<<'\n';
  	dfs<< nyr   		<<'\n';
  	dfs<< sage  		<<'\n';
  	dfs<< nage  		<<'\n';
  	dfs<< ngear 		<<'\n';
 
  	dfs<<"#Allocation"	<<'\n';
  	dfs<< dAllocation 	<<'\n';
  	
  	dfs<<"#Age-schedule and population parameters"<<'\n';
  	dfs<< d_linf  			<<'\n';
  	dfs<< d_vonbk  			<<'\n';
  	dfs<< d_to  			<<'\n';
  	dfs<< d_a  				<<'\n';
  	dfs<< d_b  				<<'\n';
  	dfs<< d_ah  			<<'\n';
  	dfs<< d_gh  			<<'\n';
  	dfs<< n_MAT				<<'\n';
	dfs<< d_maturityVector <<'\n';

  	dfs<<"#Observed catch data"<<'\n';
  	dfs<< nCtNobs 		<<'\n';
  	dfs<< dCatchData    <<'\n';

  	dfs<<"#Abundance indices"	<<'\n';
  	dfs<< nItNobs 					<<'\n';
  	dfs<< n_it_nobs 				<<'\n';
  	dfs<< n_survey_type 			<<'\n';
  	dfs<< d3_survey_data 			<<'\n';

  	dfs<<"#Age composition"		<<'\n';
  	dfs<< nAgears				<<'\n';
  	dfs<< n_A_nobs				<<'\n';
  	dfs<< n_A_sage				<<'\n';
  	dfs<< n_A_nage				<<'\n';
  	dfs<< inp_nscaler 			<<'\n';
  	dfs<< d3_A					<<'\n';

  	dfs<<"#Empirical weight-at-age data"	<<'\n';
  	dfs<< nWtTab 				<<'\n';
  	dfs<< nWtNobs				<<'\n';
	dfs<< d3_inp_wt_avg			<<'\n'; // not sure if this shoud be d3_inp_wt_avg, and how this would affect simDatfile 

	dfs<<"#EOF"	<<'\n';
	dfs<< 999	<<'\n';
	
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
FUNCTION dvector ifdSelex(const dvector& va, const dvector& ba, const double& mpow)
  {

  	dvector pa(sage,nage);

  	pa = (elem_prod(va,pow(ba,mpow)));
  	pa = pa/sum(pa);
  	pa = exp( log(pa) - log(mean(pa)) );
  	return (pa);
  }

REPORT_SECTION

	if(verbose){
    LOG<<"Start of Report Section...\n";
  }
	report<<"ObjectiveFunction\n"<<objfun<<'\n';
  report<<"FuncEvals\n"<<nf<<'\n';
  report<<"NumParams\n"<<npar<<'\n';
  report<<"MaxGrad\n"<<objective_function_value::gmax<<'\n';
  report<<"ExitCode\n"<<iexit<<'\n';
  report<<"HangCode\n"<<ihang<<'\n';
  report<<"Runtime\n"<<(long(difftime(finish,start))%3600)%60<<'\n';

	report<<DataFile<<'\n';
	report<<ControlFile<<'\n';
	report<<ProjectFileControl<<'\n';
	REPORT(objfun);
	if(!delaydiff) REPORT(nlvec);
	if(delaydiff) REPORT(nlvec_dd);
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
	report<<"rho\n"<<theta(6)<<'\n';
	report<<"vartheta\n"<<theta(7)<<'\n';
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
	//if(!delaydiff) REPORT(wa);
  	//if(!delaydiff) REPORT(fec);
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
		report<<"Neff"<<'\n';
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
					report<<sum(O(j))<<"\t"<<effectiveN<<'\n';
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

	report<<"d3_wt_avg"<<'\n';
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
			report<<d3_wt_avg(ig)(i)<<'\n';
		}
	
	}

	// |---------------------------------------------------------------------------------|
	// | SELECTIVITIES (4darray)
	// |---------------------------------------------------------------------------------|
	// |
  report<<"sel_par"<<'\n';
	for(k=1;k<=ngear;k++)
	{
	  for (j=1;j<=jsel_npar(k);j++)
	  {
	  	report<<k<<"\t"<<j<<"\t"<<exp(sel_par(k)(j))<<'\n';
	  }
	}

	report<<"log_sel"<<'\n';
	for(k=1;k<=ngear;k++)
	{
		for(int ig=1;ig<=n_ags;ig++)
		{
			for(i=syr;i<=nyr;i++)
			{
				report<<k<<"\t"<<ig<<"\t"<<i<<"\t"<<log_sel(k)(ig)(i)<<'\n';	
			}
		}
	}

	// |---------------------------------------------------------------------------------|
	// | MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	// REPORT(ft);
	report<<"ft"<<'\n';
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		report<<ft(ig)<<'\n';
	}
	report<<"ut"<<'\n';
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		report<<1.0-exp(-ft(ig))<<'\n';
	}
	if(!delaydiff){
		REPORT(M);
		REPORT(F);
		REPORT(Z);
	}else{
		REPORT(M_dd);
		REPORT(F_dd);
		REPORT(Z_dd);
		REPORT(numbers);
		REPORT(biomass);
		REPORT(annual_mean_wt);
		REPORT(wbar)
	}
	

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

  report<<"vbt"<<'\n';
		for(k=1;k<=ngear;k++){
			for(int ig=1;ig<=ngroup;ig++){
				for(i=syr;i<=nyr+1;i++){
					report<<k<<"\t"<<ig<<"\t"<<i<<"\t"<<vbt(ig)(k)(i)<<'\n';
				}
			}
	 }

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
		LOG<<"Calculating MSY-based reference points\n";
		calcReferencePoints();
		LOG<<"Finished calcReferencePoints\n";
		//exit(1);
		REPORT(bo);
		//REPORT(fmsy);
		//REPORT(msy);
		//REPORT(bmsy);
		//AQUI
		REPORT(fmsy);
		REPORT(msy);
		REPORT(bmsy);
		// REPORT(Umsy);
    LOG<<"Running Projections\n";
    //RF RE-INSTATED PROJECTION_MODEL :: ONLY IMPLEMENTED FOR AGS=1 AND FOR GEAR 1 (FISHERY)
    if(n_ags==1) {
		  int ii;
		  for(ii=1;ii<=n_tac;ii++){
        //LOG<<ii<<" "<<tac(ii)<<'\n';
		   	if (!delaydiff) projection_model(tac(ii));
		   	//if(delaydiff) 	projection_model_dd(tac(ii));
		 	}
		 }
		 if(n_ags>1){
       if(nf==1) LOG<<"************Projections not yet implemented for number of areas/groups > 1************\n\n";
		 }
     LOG<<" ______________________________ \n";
     LOG<<"|    END OF REPORT SECTION     |\n";
	   LOG<<"|______________________________|\n";
	}

	// |---------------------------------------------------------------------------------|
	// | OUTPUT FOR OPERATING MODEL
	// |---------------------------------------------------------------------------------|
	// | Move to final section?
	if( last_phase() )
	{
		ofstream ofs("iSCAM.res");

		ofs<<"# Bo\n"<<bo<<'\n';
		ofs<<"# Fmsy\n"<<fmsy<<'\n';
		ofs<<"# MSY\n"<<msy<<'\n';
		ofs<<"# Bmsy\n"<<bmsy<<'\n';
		ofs<<"# Sbt\n";
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<sbt(g)(nyr+1)<<"\t";
		}
		ofs<<'\n';

		// projected biomass
		// The total biomass for each stock
		ofs<<"# Total biomass\n";
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<bt(g)(nyr+1)<<"\t";
		}
		ofs<<'\n';

		ofs<<"# Numbers-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<N(ig)(nyr+1)<<'\n';
		}

		ofs<<"# Weight-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<d3_wt_avg(ig)(nyr+1)<<'\n';
		}		

		ofs<<"# Natural mortality-at-age\n";
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			ofs<<M(ig)(nyr)<<'\n';
		}		


		// 4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
		ofs<<"# log_selectivity\n";
		for(int k = 1; k <= ngear; k++ )	
		{
			for(int ig = 1; ig <= n_ags; ig++ )
			{
				ofs<<log_sel(k)(ig)(nyr)<<'\n';
			}
		}
	}


// 	/*
// 	Stock status info
// 	Bstatus = sbt/bmsy;
// 	Fstatus = ft/fmsy; If fmsy > 0 
// 	*/
// 	if(bmsy>0)
// 	{
// 		dvector Bstatus=value(sbt/bmsy);
// 		REPORT(Bstatus);
// 	}
	
// 	dmatrix Fstatus(1,ngear,syr,nyr);
// 	Fstatus.initialize();
// 	for(k = 1; k <= nfleet; k++)
// 	{
// 		if(fmsy(k) >0 )
// 		{
// 			j    = nFleetIndex(k);
// 			Fstatus(j) = value(ft(j)/fmsy(k));
// 		}
// 	}
// 	REPORT(Fstatus);
	
// 	//Parameter controls
// 	dmatrix ctrl=theta_control;
// 	REPORT(ctrl);
	
	
// 	if(last_phase()) decision_table();
	
	
// 	dvector rt3(1,3);
// 	if(last_phase())
// 	{
// 		dvector rt3 = age3_recruitment(value(column(N,3)),d3_wt_avg(nyr+1,3),value(M_tot(nyr,3)));
// 		REPORT(rt3);
// 	}
	
// 	//dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),d3_wt_avg(nyr+1)));
// 	dvector future_bt = value(elem_prod(N(nyr+1)*exp(-m_bar),d3_wt_avg(nyr+1)));
// 	REPORT(future_bt);
// 	double future_bt4 = sum(future_bt(4,nage));
// 	REPORT(future_bt4);
	

	if(verbose){
    LOG<<"END of Report Section...\n";;
  }

//FUNCTION decision_table
//   {
 	/*
 	This function takes a vector of projected catches and computes the following
 	Reference points: Bmsy, Bo, Fmsy, Umsy.
	
 	Biomass Metrics for the decision table:
 	1) P(SB_{t+1} < SB_{t})
 	2) P(SB_{t+1} < 0.25 B_{0})
 	3) P(SB_{t+1} < 0.75 B_{0})
 	4) P(SB_{t+1} < 0.40 B_{MSY})
 	5) P(SB_{t+1} < 0.80 B_{MSY})
	
 	Harvest Metrics for the decision table:
 	1) P(U_{t+1} > Target harvest rate)
 	2) P(U_{t+1} > 1/2 Fmsy)
 	3) P(U_{t+1} > 2/3 Fmsy)
 	4) P(tac/3+  > 20%)
	
 	Key to the harvest metric is the definition of Umsy and dAllocation to fleets.
	
 	Pseudocode:
 		1) Calculate reference points (Fmsy, Bmsy)
 		2) Loop over vector of proposed catches
 		3) Evaluate biomass metrics for each posterior sample
 		4) Evaluate harvest metrics for each posterior sample
	
 	*/
 	//int i;
 	//for(i=1;i<=n_tac;i++)
 	//{
 	//	LOG<<i<<" "<<tac<<'\n';
 	//	projection_model(tac(i));
 	//}
 	// LOG<<"Ok to here"<<'\n';			       a
     //}

FUNCTION mcmc_output
  int iter;
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
    for(int fleet=1;fleet<=nfleet;fleet++){
      ofs<<","<<"umsy"<<fleet;
    }
    for(int i=1;i<=nItNobs;i++){
      ofs<<","<<"q"<<i;
    }
    for(int group=1;group<=ngroup;group++){
      ofs<<","<<"SSB"<<group;
    }
    for(k=1;k<=ngear;k++){
      for (j=1;j<=jsel_npar(k);j++){
        ofs<<","<<"sel_g"<<k;
        ofs<<","<<"sel_sd"<<k;
      }
    }
    ofs<<","<<"f";
    ofs<<'\n';

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
    of1<<'\n';

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
    of2<<'\n';

    ofstream of3("iscam_ft_mcmc.csv");
    iter = 1;
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
    of3<<'\n';

    ofstream of4("iscam_rdev_mcmc.csv");
    iter = 1;
    for(int ag=1;ag<=n_ag;ag++){
      for(int yr=syr;yr<=nyr;yr++){
        if(iter == 1){
          of4<<"rdev"<<ag<<"_"<<yr;
        }else{
          of4<<",rdev"<<ag<<"_"<<yr;
        }
        iter++;
      }
    }
    of4<<'\n';

    ofstream of5("iscam_vbt_mcmc.csv");
    iter = 1;
    for(int ag=1;ag<=ngroup;ag++){
      for(int gear=1;gear<=ngear;gear++){
        for(int yr=syr;yr<=nyr+1;yr++){
          if(iter == 1){
            of5<<"vbt"<<ag<<"_gear"<<gear<<"_"<<yr;
          }else{
            of5<<",vbt"<<ag<<"_gear"<<gear<<"_"<<yr;
          }
          iter++;
        }
      }
    }
    of5<<'\n';

    ofstream of6("iscam_ut_mcmc.csv");
    iter = 1;
    for(int ag=1;ag<=n_ags;ag++){
      for(int gear=1;gear<=ngear;gear++){
        for(int yr=syr;yr<=nyr;yr++){
          if(iter == 1){
            of6<<"ut"<<ag<<"_gear"<<gear<<"_"<<yr;
          }else{
            of6<<",ut"<<ag<<"_gear"<<gear<<"_"<<yr;
          }
          iter++;
        }
      }
    }
    of6<<'\n';
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
  for(int fleet=1;fleet<=nfleet;fleet++){
    ofs<<","<<1.0-exp(-fmsy(fleet));
  }
  for(int it=1;it<=nItNobs;it++){
    ofs<<","<<q(it);
  }
  for(int group=1;group<=ngroup;group++){
    ofs<<","<<sbt(group)(nyr);
  }
  for(k=1;k<=ngear;k++){
    for (j=1;j<=jsel_npar(k);j++){
      ofs<<","<<exp(sel_par(k)(j)(1));
      ofs<<","<<exp(sel_par(k)(j)(2));
    }
  }

  ofs<<","<<objfun;
  ofs<<'\n';

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
  of1<<'\n';

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
  of2<<'\n';

  // output fishing mortality
  ofstream of3("iscam_ft_mcmc.csv",ios::app);
  iter = 1;
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
  of3<<'\n';

  // output recruitment deviations
  // This is what the declaration of log_dev_recs looks like:
  // init_bounded_matrix log_rec_devs(1,n_ag,syr,nyr,-15.,15.,2);
  ofstream of4("iscam_rdev_mcmc.csv",ios::app);
  iter = 1;
  for(int ag=1;ag<=n_ag;ag++){
    for(int yr=syr;yr<=nyr;yr++){
      if(iter == 1){
        of4<<log_rec_devs(ag)(yr);
      }else{
        of4<<","<<log_rec_devs(ag)(yr);
      }
      iter++;
    }
  }
  of4<<'\n';

  // output vulnerable biomass to all gears //Added by RF March 19 2015
  ofstream of5("iscam_vbt_mcmc.csv",ios::app);
  iter = 1;
  for(int ag=1;ag<=ngroup;ag++){
    for(int gear=1;gear<=ngear;gear++){
      for(int yr=syr;yr<=nyr;yr++){
        if(iter == 1){
          of5<<vbt(ag)(gear)(yr);
        }else{
          of5<<","<<vbt(ag)(gear)(yr);
        }
        iter++;
      }
    }
  }
  of5<<'\n';

  // output fishing mortality as U (1-e^-F)
  ofstream of6("iscam_ut_mcmc.csv",ios::app);
  iter = 1;
  for(int ag=1;ag<=n_ags;ag++){
    for(int gear=1;gear<=ngear;gear++){
      for(int yr=syr;yr<=nyr;yr++){
        if(iter == 1){
          of6<<1.0-exp(-ft(ag)(gear)(yr));
        }else{
          of6<<","<<1.0-exp(-ft(ag)(gear)(yr));
        }
        iter++;
      }
    }
  }
  of6<<'\n';

  ofs.flush();
  of1.flush();
  of2.flush();
  of3.flush();
  of4.flush();
  of5.flush();
  of6.flush();

 //RF:: March 17 2015. RF re-instated projection_model for Arrowtooth Flounder assessment. NOT IMPLEMENTED FOR MULTIPLE AREA/GROUPS
 // CW: Took this out while testing  the mltiple area delaydiff
 //if(n_ags==1) {
 //  int ii;
 //  for(ii=1;ii<=n_tac;ii++){
 //    //LOG<<ii<<" "<<tac(ii)<<'\n';
 //    projection_model(tac(ii));
 //  }
 //}
 //if(n_ags>1){
 //  if(nf==1) LOG<<"************Decision Table not yet implemented for number of areas/groups > 1************\n\n";
 //}

// FUNCTION dvector age3_recruitment(const dvector& rt, const double& wt,const double& M)
//   {
// /*
// This routine returns the poor average and good age-3 recruits
// that is used in constructing the decision table for the pacific
// herring fisheries.

// -1) sort the rt vector from small to large
// -2) find the 33rd and 66th percentiles
// -3) take the averge of the 0-33, 34-66, 67-100
// -4) multiply by the average weight
// -5) return the age-3 recruitment biomass
// */

// dvector s_rt = sort(rt);
// dvector rbar(1,3);

// double idx = floor((nyr-syr+1.)/3.);
// int ix1 = syr+int(idx);
// int ix2 = syr+int(2.*idx);
// rbar(1) = mean(s_rt(syr,ix1));
// rbar(2) = mean(s_rt(ix1+1,ix2));
// rbar(3) = mean(s_rt(ix2+1,nyr));
// rbar = rbar*wt*exp(-M);
// //LOG<<rbar<<'\n';
// return(rbar);
//   }

 //RF re-instated this code (with several updates to match current version, and new code for projection output files for 2014 Arrowtooth Flounder) March 17 2015
 // !!! NOT IMPLEMENTED FOR n_ags > 1!!!
FUNCTION void projection_model(const double& tac);
  /*
  This routine conducts population projections based on
  the estimated values of theta.  Note that all variables
  in this routine are data type variables.

  Arguments:
  tac is the total allowable catch that must be allocated
  to each gear type based on dAllocation(k)

  theta(1) = log_ro
  theta(2) = h
  theta(3) = log_m
  theta(4) = log_avgrec
  theta(5) = log_recinit
  theta(6) = rho
  theta(7) = vartheta

  ** NOTES **
  * Projections are based on average natural mortality and fecundity.
  * Selectivity is based on selectivity in terminal year.
  * Average weight-at-age is based on mean weight in the last 5 years.
  */
  static int runNo=0;
  runNo ++;
  int i,k;
  int pyr = nyr+1;
  BaranovCatchEquation cBaranov;
  // | (2) : Average weight and mature spawning biomass for reference years  (copied from calcReferencePoints() but only implemented for ig=1)
  // |     : dWt_bar(1,n_ags,sage,nage)
  dvector fa_bar(sage,nage);
  dvector  M_bar(sage,nage);
  fa_bar = elem_prod(dWt_bar(1),ma(1));
  M_bar  = colsum(value(M(1).sub(pf_cntrl(3),pf_cntrl(4))));
  M_bar /= pf_cntrl(4)-pf_cntrl(3)+1;
  // --derive stock re4cruitment parameters
  // --survivorship of spawning biomass
  dvector lx(sage,nage);
  double  tau = value(sqrt(1.-rho)*varphi);
  //double m_rho = d_iscamCntrl(13);
  lx(sage)     = 1.;
  for(i=sage+1; i<=nage; i++){
   lx(i) = lx(i-1)*mfexp(-M_bar(i-1));
   if(i==nage)  lx(i) /= 1.0 - mfexp( -M_bar(i));
  }
  double phib = lx*fa_bar;   //average fecundity is calculated for all area/groups but projections are currently only implemented for n_ags=1
  double so = value(kappa(1)/phib);
  double bo  = value(ro(1)*phib);
  double beta = 1;
  switch(int(d_iscamCntrl(2))){
   case 1:  // Beverton-Holt
    beta = value((kappa(1)-1.)/bo);
    break;
   case 2:  // Ricker
    beta = value(log(kappa(1)/bo));
    break;
  }
  /* Fill arrays with historical values */
  dvector p_sbt(syr,pyr+1);
  dvector  p_ct(1,ngear);
  dmatrix  p_ft(nyr,pyr+1,1,ngear);
  dmatrix   p_N(syr,pyr+2,sage,nage);
  dmatrix   p_Z(syr,pyr+1,sage,nage);
  p_N.initialize();
  p_sbt.initialize();
  p_Z.initialize();
  p_ct.initialize();
  p_ft.initialize();

  //The main model already does a projection to nyr+1
  //but want to draw an average recruitment for projection rather than highly uncertain estimate
  for(i = syr; i<=nyr; i++){
   p_N(i) = value(N(1)(i));
   p_sbt(i) =  value(sbt(1)(i));
   p_Z(i) =  value(Z(1)(i));
  }

  /* Selectivity and dAllocation to gears */
  dmatrix va_bar(1,ngear,sage,nage);
  for(k=1;k<=ngear;k++){
   p_ct(k)   = dAllocation(k)*tac;
   va_bar(k) = exp(value(log_sel(k)(1)(nyr)));
  }

  /* Simulate population into the future under constant tac policy. */
  for(i = nyr-1; i<=pyr+1; i++){
    //ft(nyr) is a function of ct(nyr) not the tac so use ft(nyr) from the main model for nyr
    if(i>nyr){
     // get_ft is defined in the Baranov.cpp file
     //average weight is calculated for all area/groups but projections are currently only implemented for n_ags=1
     p_ft(i) = cBaranov.getFishingMortality(p_ct,M_bar, va_bar, p_N(i),dWt_bar(1));  
     // calculate total mortality in future years
     p_Z(i) = M_bar;
     for(k=1;k<=ngear;k++){
      p_Z(i)+=p_ft(i,k)*va_bar(k);
     }
    }

    //Overwrite sbt(nyr) so that it does not include estimated Rt(nyr), which is highly uncertain
    //This will only be different from sbt(nyr) in the main model if recruits contribute to the spawning population, which is rare
    //d_iscamCntrl(13) is defined as: fraction of total mortality that takes place prior to spawning
    if(i>=nyr){
      p_sbt(i) = elem_prod(p_N(i),exp(-p_Z(i)*d_iscamCntrl(13))) *fa_bar;
    }

    // sage recruits with random deviate xx
    // note the random number seed is repeated for each tac level.
    //NOTE that this treatment of rec devs is different from historical model
    double  xx = randn(nf+i)*tau;
    if(i>=syr+sage-1){
      double rt = 1;
      double et = p_sbt(i-sage+1);  //lagged spawning biomass  (+1 because we want recruits for year i+1)
      if(d_iscamCntrl(2)==1){      // Beverton-Holt model
        rt=(so*et/(1.+beta*et));
      }
      if(d_iscamCntrl(2)==2){      // Ricker model
        rt=(so*et*exp(-beta*et));
      }
      p_N(i+1,sage)=rt*exp(xx-0.5*tau*tau);  //Next year's recruits
    }
    /* Update numbers at age in future years*/
    //Next year's numbers
    p_N(i+1)(sage+1,nage) =++ elem_prod(p_N(i)(sage,nage-1),exp(-p_Z(i)(sage,nage-1)));
    p_N(i+1,nage)        +=   p_N(i,nage)*exp(-p_Z(i,nage));

    //Predicted catch for checking calculations    (RF tested this March 17, 2015)
    //if(i > nyr){
    //        LOG<<p_ft<<'\n';
    // for(k=1;k<=ngear;k++)
    // {
    //  dvector ba = elem_prod(p_N(i),dWt_bar(1));
    //  double ctest;
    //  ctest = sum(elem_div(elem_prod(elem_prod(ba,p_ft(i,k)*va_bar(k)),1.-exp(-p_Z(i))),p_Z(i)));
    //  LOG<<"gear = "<<k<<" tac = "<<tac<<"\t ct = "<<ctest<<'\n';
    // }
    //}
  } //end year loop
  /*
  Write output to projection file for constructing decision tables.
  For BC Arrowtooth Flounder 2014 assessment (Forrest, Grandin, Pacific Biological Station)
  yr  = 2014; pyr = 2015; pyr+1 = 2016
  the object p_ft only has one year of ft projections (pyr), but has ngears.
  FOR THE ARROWTOOTH ASSESSMENT ONLY WRITE OUT THE FIRST GEAR
  */

  if(mceval_phase()){
   if(nf==1 && runNo==1){
    LOG<<"Running MCMC projections\n";
    ofstream ofsmcmc("iscammcmc_proj_Gear1.csv");
    write_proj_headers(ofsmcmc, syr, nyr);
    ofsmcmc.flush();
   }
   ofstream ofsmcmc("iscammcmc_proj_Gear1.csv", ios::app);
   write_proj_output(ofsmcmc, syr, nyr, tac, pyr, p_sbt, p_ft, ft(1), bo, fmsy, bmsy);
   ofsmcmc.flush();
  }else{
//   if(runNo==1){
//    LOG<<"Running MPD projections"<<'\n';
//    ofstream ofsmpd("iscammpd_proj_Gear1.csv");
//    write_proj_output(ofsmpd);
//    ofsmpd.flush();
//   }
//   ofstream ofsmpd("iscammpd_proj_Gear1.csv", ios::app);
//   write_proj_output(ofsmpd, tac, pyr, p_sbt, p_ft);
//   ofsmpd.flush();
  }
  if(!mceval_phase()){
   LOG<<"Finished projection model for TAC = "<<tac<<'\n';
  }



FUNCTION void projection_model_dd(const double& tac);	
  {
	/*
	This routine conducts population projections based on 
	the estimated values of theta.  Note that all variables
	in this routine are data type variables.
	
	Arguments:
	tac is the total allowable catch that must be allocated 
	to each gear type based on allocation(k)
	
	theta(1) = log_ro
	theta(2) = h
	theta(3) = log_m
	theta(4) = log_avgrec
	theta(5) = log_recinit
	theta(6) = rho
	theta(7) = vartheta
	
	** NOTES **
	* Projections are based on estimated constant natural mortality 
	
	*/
	static int runNo=0;
	runNo ++;
	int i;
	int pyr = nyr+2;	//projection year. 

	 BaranovCatchEquation cBaranov;
	     
	//get parameters - convert to data objects
	//double pbo   = value(bo(1));
	//double pso = value (so(1));
	//double pbeta =value(beta(1));
		  
	dvector p_bt(syr,pyr);
	dvector p_ft(syr,pyr);
	dvector p_N(syr,pyr);
	dvector p_S(syr,pyr);
	dvector p_rt(syr+sage,pyr);//dmatrix p_Ft(1,ngear,syr,pyr);
	p_bt.initialize();
	p_ft.initialize();
	p_N.initialize();
	p_S.initialize();
	p_rt.initialize();
	
	p_ft(syr,nyr) = value(ft(1)(1)(syr,nyr));
	p_N(syr,nyr) = value(numbers(1)(syr,nyr));
	p_bt(syr,nyr)   = value(biomass(1)(syr,nyr)); //sbt and vul biomass all the same for delay diff
	p_S(syr,nyr)   = value(surv(1)(syr,nyr));
	p_rt(syr+sage,nyr)   = value(rt(1)(syr+sage,nyr));
		
	//control points    - these are "historical" control points based on biomass and F reconstruction
	// Question!! I am not sure if these numbers should be fixed, maybe read from the pfc file?


	int nshort=pf_cntrl(7)-syr+1;
	int nlong=pf_cntrl(8)-syr+1;
	double meanfshort;	  // average F between 1956 and 2004
	double meanflong;	    // average F between 1956 and 2012
	 double meanbshort;	 // average B between 1956 and 2004
	double meanblong;	  // average B between 1956 and 2012
	double minb;	  // biomass in 1971 for 5CD or 1985 for 5AB

	dvector hist_ftshort(syr,pf_cntrl(7));
	dvector hist_ftlong(syr,pf_cntrl(8));
	dvector hist_btshort(syr,pf_cntrl(7));
	dvector hist_btlong(syr,pf_cntrl(8));
    hist_ftshort.initialize();  hist_ftlong.initialize();
	hist_btshort.initialize();  hist_btlong.initialize();

		

	hist_ftshort=value(ft(1)(1)(syr,pf_cntrl(7)));
	hist_btshort=value(biomass(1)(syr,pf_cntrl(7)));
	
	if(nyr>=pf_cntrl(8)){

		hist_ftlong=value(ft(1)(1)(syr,pf_cntrl(8)));
		hist_btlong=value(biomass(1)(syr,pf_cntrl(8)));
	}


	meanfshort=sum(hist_ftshort)/nshort;
	if(nyr>=pf_cntrl(8)) meanflong=sum(hist_ftlong)/nlong;
	meanbshort=sum(hist_btshort)/nshort;
	if(nyr>=pf_cntrl(8)) meanblong=sum(hist_btlong)/nlong;
	

	// Question CW - Where does 1985 comes from, it it the minimum biomass observed or is it a set number?
	minb=hist_btshort(1985);

	/* Simulate population into the future under constant tac policy. */
	
	for(i = nyr+1; i<=pyr; i++)
	{
		//recruits
		//double p_tau = value(sqrt(1-rho)/varphi);
		double p_tau = value(tau(1)); 

		//question CW What's nf??
		double xx = randn(nf+i)*p_tau;
			
		

		//double rt;
		double et=p_bt(i-kage(1)); //delay diff

		if(d_iscamCntrl(2)==1)p_rt(i)=value((so(1)*et/(1.+beta(1)*et))*exp(xx-0.5*p_tau*p_tau));
		if(d_iscamCntrl(2)==2)p_rt(i)=value((so(1)*et*exp(-beta(1)*et))*exp(xx-0.5*p_tau*p_tau));
				
		//numbers and biomass
		//Update biomass and numbers	
		p_bt(i) =(p_S(i-1)*(rho_g(1)*p_bt(i-1)+alpha_g(1)*p_N(i-1))+wk(1)*p_rt(i));
		p_N(i)=p_S(i-1)*p_N(i-1)+p_rt(i);
               	//LOG<<i<<"  "<<p_rt(i)<<"  "<<p_bt(i)<<"  "<<p_N(i)<<'\n';

		//get_ft is defined in the Baranov.cpp file
		p_ft(i) = cBaranov.get_ftdd(tac,value(M_dd(1)(syr)),p_bt(i));	    //hardwiring the catch to gear 1 for this assessment       m_bar is same as constant M

		/*
		//test get_ftdd with Baranov equation
		double testf =p_ft(i);
		double testc = p_bt(i)*(1-mfexp(-value(m)-testf))*(testf/(value(m) + testf));
		LOG<<i<<"  "<<testf<<"  "<<tac<<"  "<<testc<<'\n'<<'\n';
		*/
		
		//Calculate mortality for next projection year
		p_S(i) = mfexp(-(value(M_dd(1)(syr))+p_ft(i)));
				
	} 
	
	  //S= Short (1956-2004) 
	   //L-Long(1956-2012)
	if(mceval_phase()){
		if(nf==1 && runNo==1)
		{
			ofstream ofsP("iscammcmc.proj");
			ofsP<<"tac" <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+1 <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<"B_"<<nyr+1 <<setw(6)     <<   "\t";		   //want probability B2015<B2014 - this will be < 1 if true
			ofsP<<"F_"<<nyr <<setw(6)     <<   "\t";
			ofsP<<"F_"<<nyr+1 <<setw(6)     <<   "\t";
			ofsP<<"F_"<<nyr+1<<"F_"<<nyr <<setw(6)     <<   "\t";		   //want probability F2014>F2013     - this will be > 1 if true
			//MSY based ref points
			ofsP<<"BMSY" <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<"BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<BMSY - this will be < 1 if true
			ofsP<<"B_"<<nyr+2<<"0.8BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.8BMSY - this will be< 1 if true
			ofsP<<"B_"<<nyr+2<<"0.4BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.4BMSY - this will be < 1 if true
			ofsP<<"FMSY" <<setw(6)     <<   "\t";
			ofsP<<"F_"<<nyr+1<<"FMSY"<<setw(6)     <<   "\t";		   //want probability F2014>F2013 - this will be > 1 if true
			//Historical ref points "short"	 1956-2004
			ofsP<<"Bmin" <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<"Bmin" <<setw(6)     <<   "\t";		   //want probability B2015<Bmin 
			ofsP<<"BAvg_S" <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<"BAvg_S" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg 
			ofsP<<"FAvg_S" <<setw(6)     <<   "\t";
			ofsP<<"F_"<<nyr+1<<"FAvg_S"<<setw(6)     <<   "\t";	
			//Historical ref points "long"	 1956-2012
			ofsP<<"BAvg_L" <<setw(6)     <<   "\t";
			ofsP<<"B_"<<nyr+2<<"BAvg_L" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg - this will be < 1 if true
			ofsP<<"FAvg_L" <<setw(6)     <<   "\t";
			ofsP<<"F_"<<nyr+1<<"FAvg_L\n";		   //want probability F2014>F2013 - this will be > 1 if true
		      
			LOG<<"Running MCMC evaluations"<<'\n';
			LOG<<"Bo when nf==1 \t"<<bo<<'\n';
		}

		ofstream ofsP("iscammcmc.proj",ios::app);
		ofsP <<tac <<setw(6)                            <<"\t"
		  << p_bt(pyr-1) <<setw(6)       <<"\t"	      
		  << p_bt(pyr) <<setw(6)       <<"\t"		 
		  << p_bt(pyr)/p_bt(pyr-1) <<setw(6)      <<"\t"	     
		 << p_ft(pyr-2) <<setw(6)      <<"\t"
		 << p_ft(pyr-1)  <<setw(6)     <<"\t"
		 << p_ft(pyr-1)/p_ft(pyr-2)  <<setw(6)     <<"\t"	 
		//MSY based ref points
		<<bmsy <<setw(6)     <<   "\t"
		<<p_bt(pyr)/bmsy <<setw(6)     <<   "\t"		 
		<<p_bt(pyr)/(0.8*bmsy) <<setw(6)     <<   "\t"		  
		<<p_bt(pyr)/(0.4*bmsy) <<setw(6)     <<   "\t"		   
		<<fmsy <<setw(6)     <<   "\t"
		<<p_ft(pyr-1)/fmsy <<setw(6)     <<   "\t"		   
		//Historical ref points "short"	 1956-2004
		<<minb <<setw(6)     <<   "\t"
		<<p_bt(pyr)/minb <<setw(6)     <<   "\t"		   
		<<meanbshort <<setw(6)     <<   "\t"
		<<p_bt(pyr)/meanbshort <<setw(6)     <<   "\t"		   
		<<meanfshort <<setw(6)     <<   "\t"
		<<p_ft(pyr-1)/meanfshort<<setw(6)     <<   "\t"		  
		 //Historical ref points "long"	 1956-2012
		<<meanblong <<setw(6)     <<   "\t"
		<<p_bt(pyr)/meanblong <<setw(6)     <<   "\t"		   
		<<meanflong <<setw(6)     <<   "\t"
		<<p_ft(pyr-1)/meanflong<<   "\t"		   	   		   
		 <<'\n';
	   }

	   //S= Short (1956-2004) 
	   //L-Long(1956-2012)
	   if(last_phase() && !mceval_phase()){
	   		if(runNo==1)
	   		{
	   			LOG<<"Running MPD projections"<<'\n';
	   			
	   			ofstream ofsP("iscammpd.proj");
	   			ofsP<<"tac" <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+1 <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2 <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2<<"B_"<<nyr+1 <<setw(6)     <<   "\t";		   //want probability B2015<B2014 - this will be < 1 if true
				ofsP<<"F_"<<nyr <<setw(6)     <<   "\t";
				ofsP<<"F_"<<nyr+1 <<setw(6)     <<   "\t";
				ofsP<<"F_" <<nyr+1<<"F_"<<nyr <<setw(6)     <<   "\t";		   //want probability F2014>F2013     - this will be > 1 if true
				//MSY based ref points
				ofsP<<"BMSY" <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2<<"BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<BMSY - this will be < 1 if true
				ofsP<<"B_"<<nyr+2<<"0.8BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.8BMSY - this will be< 1 if true
				ofsP<<"B_"<<nyr+2<<"0.4BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.4BMSY - this will be < 1 if true
				ofsP<<"FMSY" <<setw(6)     <<   "\t";
				ofsP<<"F_"<<nyr+1<<"FMSY"<<setw(6)     <<   "\t";		   //want probability F2014>F2013 - this will be > 1 if true
				//Historical ref points "short"	 1956-2004
				ofsP<<"Bmin" <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2<<"Bmin" <<setw(6)     <<   "\t";		   //want probability B2015<Bmin 
				ofsP<<"BAvg_S" <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2<<"BAvg_S" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg  
				ofsP<<"FAvg_S" <<setw(6)     <<   "\t";
				ofsP<<"F_"<<nyr+1<<"FAvg_S"<<setw(6)     <<   "\t";
				//Historical ref points "long"	 1956-2012
				ofsP<<"BAvg_L" <<setw(6)     <<   "\t";
				ofsP<<"B_"<<nyr+2<<"BAvg_L" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg - this will be < 1 if true
				ofsP<<"FAvg_L" <<setw(6)     <<   "\t";
				ofsP<<"F_"<<nyr+ 1<< "FAvg_L\n";		   //want probability F2014>F2013 - this will be > 1 if true
				
	   		}
	   
	   		LOG<<"tac = "<<tac<<'\n';
	   		ofstream ofsP("iscammpd.proj",ios::app);
	   		ofsP 
	   		  <<tac <<setw(6)                            <<"\t"
			  << p_bt(pyr-1) <<setw(6)       <<"\t"	      
			  << p_bt(pyr) <<setw(6)       <<"\t"		 
			  << p_bt(pyr)/p_bt(pyr-1) <<setw(6)      <<"\t"	     
		  	 << p_ft(pyr-2) <<setw(6)      <<"\t"
			 << p_ft(pyr-1)  <<setw(6)     <<"\t"
		  	 << p_ft(pyr-1)/p_ft(pyr-2)  <<setw(6)     <<"\t"	 
			//MSY based ref points
			<<bmsy <<setw(6)     <<   "\t"
			<<p_bt(pyr)/bmsy <<setw(6)     <<   "\t"		 
			<<p_bt(pyr)/(0.8*bmsy) <<setw(6)     <<   "\t"		  
			<<p_bt(pyr)/(0.4*bmsy) <<setw(6)     <<   "\t"		   
			<<fmsy <<setw(6)     <<   "\t"
			<<p_ft(pyr-1)/fmsy <<setw(6)     <<   "\t"		   
			//Historical ref points "short"	 1956-2004
			<<minb <<setw(6)     <<   "\t"
			<<p_bt(pyr)/minb <<setw(6)     <<   "\t"		   
			<<meanbshort <<setw(6)     <<   "\t"
			<<p_bt(pyr)/meanbshort <<setw(6)     <<   "\t"		   
			<<meanfshort <<setw(6)     <<   "\t"
			<<p_ft(pyr-1)/meanfshort<<setw(6)     <<   "\t"		  
			 //Historical ref points "long"	 1956-2012
			<<meanblong <<setw(6)     <<   "\t"
			<<p_bt(pyr)/meanblong <<setw(6)     <<   "\t"		   
			<<meanflong <<setw(6)     <<   "\t"
			<<p_ft(pyr-1)/meanflong<<   "\t"		   	   		   
			<<'\n';
	   }
	}








//end of projection model dd


FUNCTION void runMSE()
	LOG<<"Start of runMSE"<<'\n';

	// STRUCT FOR MODEL VARIABLES
	// ModelVariables s_mv;
	// s_mv.log_ro    = value( theta(1) );
	// s_mv.steepness = value( theta(2) );
	// s_mv.m         = value( theta(3) );
	// s_mv.log_rbar  = value( theta(4) );
	// s_mv.log_rinit = value( theta(5) );
	// s_mv.rho       = value( theta(6) );
	// s_mv.varphi    = value( theta(7) );

	// // Selectivity parameters
	// d3_array log_sel_par(1,ngear,1,jsel_npar,1,isel_npar);
	// d4_array d4_log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
	// for(int k = 1; k <= ngear; k++ )
	// {
	// 	log_sel_par(k) = value(sel_par(k));
	// 	d4_log_sel(k)  = value(log_sel(k));
	// }
	
	// s_mv.d3_log_sel_par = &log_sel_par;
	// s_mv.d4_logSel      = &d4_log_sel;

	// d3_array d3_M(1,n_ags,syr,nyr,sage,nage);
	// d3_array d3_F(1,n_ags,syr,nyr,sage,nage);
	// for(int ig = 1; ig <= n_ags; ig++ )
	// {
	// 	d3_M(ig) = value(M(ig));
	// 	d3_F(ig) = value(F(ig));
	// }

	// s_mv.d3_M = &d3_M;
	// s_mv.d3_F = &d3_F;
	// s_mv.log_rec_devs = value(log_rec_devs);
	// s_mv.init_log_rec_devs = value(init_log_rec_devs);

	// s_mv.q = value(q);
	// s_mv.sbt = value(sbt);
	// d3_array tmp_ft=value(ft);
	// s_mv.d3_ft = &tmp_ft;

	// s_mv.sbo = value(sbo);
	// s_mv.so = value(so);


	// // |-----------------------------------|
	// // | Instantiate Operating Model Class |
	// // |-----------------------------------|
	// OperatingModel om(s_mv,argc,argv);
	// om.runScenario(rseed);

	// LOG<<"DONE\n";

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << "\n";

	#undef TINY
	#define TINY 1.e-08

	#undef NA
	#define NA -99.0

	#include <time.h>
	#include <string.h>
  #include "../../include/baranov.h"
  #include "../../include/LogisticNormal.h"
  #include "../../include/LogisticStudentT.h"
	#include "../../include/msy.h"
  #include "../../include/msy.hpp"
	#include "../../include/multinomial.h"
  #include "../../include/utilities.h"
  #include "../../include/Logger.h"

	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	bool mcmcEvalPhase = 0;
	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;
	// Variables to store results from DIC calculations.
	double dicNoPar = 0;
	double dicValue = 0;

 //Extra test functions by RF to test ref points
 //Called by run_FRP() in calcReferencePoints
FUNCTION void slow_msy(dvector& ftest, dvector& ye, dvector& be, double& msy, double& fmsy, double& bmsy )
	//THIS CODE VERIFIES THAT THE EQM CODE IS RETURNING CORRECT REF POINTS
	int i;
	int j;
	int k;
	int t;
	int NF=size_count(ftest);
	int Nyr=100; //number of years to run out the model
	ye.initialize();
	be.initialize();
	
	double sa;
	dvector za(sage,nage); za.initialize();
	dvector saf(sage,nage); saf.initialize();
	dvector lx(sage,nage); lx.initialize();
	dvector vd(sage,nage); vd.initialize();
	dvector avg_wt(sage,nage); avg_wt.initialize(); 
	dvector avg_fec(sage,nage);
	dvector  M_bar(sage,nage);

        avg_wt = dWt_bar(1);
	avg_fec = elem_prod(dWt_bar(1),ma(1));
	vd = exp(value(log_sel(1)(1)(nyr)));
	
	double Ro=value(ro(1));
	//double CR=value(kappa(1));
	double Mbar;
	M_bar  = colsum(value(M(1).sub(pf_cntrl(3),pf_cntrl(4)))); //mean across years
	M_bar /= pf_cntrl(4)-pf_cntrl(3)+1;
	Mbar = mean(M_bar); //mean across ages
	
	dmatrix Nn(1,Nyr+1,sage,nage);			//Numbers at age
	dmatrix Ff(1,Nyr+1,sage,nage);			//Age-specific fishing mortality
	dmatrix Zz(1,Nyr+1,sage,nage);
	dvector Ss(sage,nage);
	dmatrix Cc(1,Nyr,sage,nage);
	dvector Ssb (1,Nyr+1);
	dvector Bb(1,Nyr+1);
	dvector Y(1,Nyr);				//predicted catch biomass
	
	//dvector finaly(1,NF);
	//dvector finalb(1,NF);
	
	//unfished
	sa=mfexp(-Mbar);
	lx(sage)=1.0;
	for(i=(sage+1); i<=nage; i++)
		lx(i)=lx(i-1)*sa;
	lx(nage)/=(1.-sa); 
		
	//Initialize model - same for all F scenarios
	for(j=sage;j<=nage;j++) Nn(1,j)=Ro*lx(j);
	Ssb(1)= sum(elem_prod(Nn(1),avg_fec));
	     
	for(k=1;k<=NF;k++){ 
		
		za=(Mbar+ftest(k)*vd);
		Ss=mfexp(-za);
		
		/*
		LOG<<"Nn "<<'\n'<<Nn(1)<<'\n';
		LOG<<"SS "<<'\n'<<Ss<<'\n';
		LOG<<"SSb "<<'\n'<<Ssb(1)<<'\n';
		LOG<<"wt "<<avg_wt<<'\n';
		*/
			
		for(t=1;t<=Nyr;t++){
			
			Nn(t+1)(sage+1,nage)=++elem_prod(Nn(t)(sage,nage-1),Ss(sage,nage-1));
			Nn(t+1,nage)+=Nn(t,nage)*Ss(nage);
			if(t==1) Nn(t+1)(sage)=Ro;
			if(t>1) Nn(t+1)(sage)=value(so(1))*Ssb(t-1)/(1.+value(beta(1))*Ssb(t-1));
			Ssb(t+1)= sum(elem_prod(Nn(t+1),avg_fec));
			
			//catch
			for(j=sage;j<=nage;j++) Cc(t,j) = ((ftest(k)*vd(j))/za(j))*(1.-exp(-(za(j))))*Nn(t,j)*avg_wt(j);
			Y(t)=sum(Cc(t));
		}//end t
		
	 ye(k)=Y(Nyr);
	 be(k)=Ssb(Nyr);
	} //end k
	
	//ye =  finaly;
	//be = finalb;
	
	//get MSY and Fmsy
	msy=max(ye);
	double mtest;
	for(k=1; k<=NF; k++)
	{
		mtest=ye(k);
		if(mtest==msy) fmsy=ftest(k);
		if(mtest==msy) bmsy=be(k);
	}
	LOG<<"Ref points from running out model\n";;
	LOG<<msy<<'\n';
	LOG<<fmsy<<'\n';
	LOG<<bmsy<<'\n';


FUNCTION void ddiff_msy(dvector& ftest, dvector& ye, dvector& be, double& msy, double& fmsy, double& bmsy )
	
	int k ;
	int NF=size_count(ftest);
	ye.initialize();
	be.initialize();
	
	
	dvariable rec_a;
	dvariable rec_b;
	
	dvariable M;

	//double M = value(m);
	rec_a=value(so(1));
	rec_b=value(beta(1));

	
	//int f,g,h,ih,gs;

	// Calculate equilibrium survivorship as function of FMSY
	for(k=1; k<=NF; k++)
	{
		double se; //survival in equilibrium
		double we; //average weight in equilibrium
				
			se = exp(-value(M_dd(1)(nyr)) - ftest(k));
			we = (se*alpha_g(1)+wk(1) *(1.-se))/(1.-rho_g(1)*se);
			
			
			be(k) = value(-1.*((-we + se*alpha_g(1) + se*rho_g(1)*we + wk(1)*rec_a*we)/(rec_b*(-we + se*alpha_g(1) + se*rho_g(1)*we)))); //Martell
			
			M = value(M_dd(1)(nyr));
			

			ye(k)   = value(be(k)*(1.0-mfexp(-ftest(k)-M))*(ftest(k)/(ftest(k)+M)));
		  	

		  	if(ye(k)<0) ye(k)=0.;
		  	if(be(k)<0) be(k)=0.;

	}
		 
	//LOG<<"ye"<<ye<<'\n';
	//LOG<<"be"<<be<<'\n';
	//LOG<<"ftest"<<ftest<<'\n';
	//exit(1);
	
	double mtest;	
	
		msy=max(ye);
			
		for(k=1; k<=NF; k++)
		{
			mtest=ye(k);
				
			if(mtest==msy){
				fmsy=ftest(k);
				bmsy=be(k);
			} 
		}
	  	
	
	LOG<<"Slow msy calcs"<<'\n';
	LOG<<"msy"<<msy<<'\n';
	LOG<<"fmsy"<<fmsy<<'\n';
	LOG<<"bmsy"<<bmsy<<'\n'; 
	
	
	
	
	




//RF's function for calling slow msy routine to test ref points
FUNCTION void run_FRP()
	//Reference points
	if(last_phase()){
	  LOG<<"\n*********Getting reference points the slow way************\n";
	  LOG<<"*******************************************\n\n";
	}
	dvector ftest(1,4001);
	ftest.fill_seqadd(0,0.01);
	ftest(1) =  21.7117; //option to put in a test value
	int Nf;
	Nf=size_count(ftest);
	double Fmsy;
	double MSY;
	double Bmsy;
	dvector Ye(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	dvector Be(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	double fmsy,msy,bmsy;
	dvector ye(1,Nf);
	dvector be(1,Nf);

	slow_msy(ftest, ye, be, msy, fmsy, bmsy);
	Fmsy=fmsy;
	MSY=msy;
	Bmsy=bmsy;
	Ye=ye;
	Be=be;

	ofstream ofsr("TEST_frp.rep");
	ofsr<<"Fmsy"<<'\n'<<Fmsy<<'\n';
	ofsr<<"MSY"<<'\n'<<MSY<<'\n';
	ofsr<<"Bmsy"<<'\n'<<Bmsy<<'\n';
	ofsr<<"ftest"<<'\n'<<ftest<<'\n';
	ofsr<<"Ye"<<'\n'<<Ye<<'\n';	
	ofsr<<"Be"<<'\n'<<Be<<'\n';	

FUNCTION void run_FRPdd()
	

	
	if(n_ags>1){
		LOG<<"MSY quantities not defined for n_ags>1"<<'\n'; 
	}else{

	dvector ftest(1,101);
	ftest.fill_seqadd(0,0.01);
	
	int Nf;
	Nf=size_count(ftest);
	

		

	for(int g=1; g<=ngroup; g++)
	{
		//dvector Ye(1,Nf); // i think this should be a matrix by group and gear
		//Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
		//dvector Be(1,Nf); // i think this should be a matrix by group
		//Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	
		dvector ye(1,Nf); // i think this should be a matrix by group and gear
		dvector be(1,Nf);

		//double fmsy;
		//double msy;
		//double bmsy;

		ddiff_msy(ftest, ye, be, msy(g,1), fmsy(g,1), bmsy(g));

		//Fmsy=fmsy;
		//MSY=msy;
		//Bmsy=bmsy;
		//Ye(g)=ye(g);
		//Be(g)=be(g);
	}


	}
	
	
	
	


	
	


FINAL_SECTION
	LOG<<"\n\nNumber of function evaluations: "<<nf<<'\n';

