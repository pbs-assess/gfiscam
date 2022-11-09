// iSCAM = integrated Statistical Catch Age Model

DATA_SECTION
	!! LOG<<"DATA_SECTION\n";
	// |---------------------------------------------------------------------------------|
	// | STRINGS FOR INPUT FILES                                                         |
	// |---------------------------------------------------------------------------------|
	// | DataFile.dat           : data to condition the assessment model on
	// | ControlFile.ctl        : controls for phases, selectivity options
	// | ProjectFileControl.pfc : used for stock projections under TAC
	// | BaseFileName           : file prefix used for all iSCAM model output
	// | ReportFileName         : file name to copy report file to.
	// |
	init_adstring DataFile;
	init_adstring ControlFile;
	init_adstring ProjectFileControl;
	init_adstring ProcedureControlFile;
	init_adstring ScenarioControlFile;
	!! BaseFileName = stripExtension(ControlFile);
	!! ReportFileName = BaseFileName + adstring(".rep");

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
	init_int n_projyr; // Number of years to project the stock
	init_int n_tac;    // Number of catch options to explore in the decision table
	init_vector tac(1,n_tac);
	init_int n_pfcntrl;
	init_vector pf_cntrl(1,n_pfcntrl);
	init_int eof_pf;
	LOC_CALCS
	  if(eof_pf != -999){
	    LOG<<"Error reading projection file.\n";
	    LOG<<"Last integer read is "<<eof_pf<<'\n';
	    LOG<<"The file should end with -999.\n Aborting!\n";
	    ad_exit(1);
	  }
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | COMMAND LINE ARGUMENTS
	// |---------------------------------------------------------------------------------|
	// | retro_yrs  : number of terminal years to remove.
	int retro_yrs;
	int testMSY;
	LOC_CALCS
	  int on, opt;
	  // command line option for retrospective analysis. "-retro retro_yrs"
	  retro_yrs = 0;
	  if((on = option_match(ad_comm::argc,ad_comm::argv, "-retro", opt)) > -1){
	    retro_yrs = atoi(ad_comm::argv[on + 1]);
	    LOG<<" ____________________________________________________ \n";
	    LOG<<"| Implementing Retrospective analysis                |\n";
	    LOG<<"|____________________________________________________|\n";
	    LOG<<"| Number of retrospective years = "<<retro_yrs<<'\n';
	  }
	END_CALCS
	// Change datafile name to be the name found on the first line in
	//  the original datafile (iscam.dat)
	!! ad_comm::change_datafile_name(DataFile);

	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS
	// |---------------------------------------------------------------------------------|
	// | area  : f
	// | group : g
	// | sex   : h
	// | year  : i
	// | age   : j
	// | gear  : k
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
	init_number propfemale;
	vector age(sage,nage);

	// |---------------------------------------------------------------------------------|
	// | LINKS TO MANAGE ARRAY INDEXING
	// |---------------------------------------------------------------------------------|
	// | n_ags    : total number of areas * groups * sex
	// | n_ag     : total number of areas * groups
	// | n_gs     : total number of groups * sex
	// | n_area   : vector of indexes for area for each sex & group combination.
	// | n_group  : vector of indexes for stock for each sex & area combination.
	// | n_sex    : vector of indexes for sex foe each area & group combination.
	// | pntr_ag  : matrix of indices for area and group.
	// | pntr_gs  : matrix of indices for group and sex.
	// | pntr_ags : d3_array of indices for area group sex.
	int n_ags;
	!! n_ags = narea * ngroup * nsex;
	int n_ag;
	!! n_ag = narea * ngroup;
	int n_gs;
	!! n_gs = ngroup * nsex;
	ivector n_area(1,n_ags);
	ivector n_group(1,n_ags);
	ivector n_sex(1,n_ags);
	imatrix pntr_ag(1,narea,1,ngroup);
	imatrix pntr_gs(1,ngroup,1,nsex);
	3darray pntr_ags(1,narea,1,ngroup,1,nsex);
	LOC_CALCS
	  age.fill_seqadd(sage,1);
	  int ig, ih, is;
	  ig = 0;
	  ih = 0;
	  is = 0;
	  for(f = 1; f <= narea; f++){
	    for(g = 1; g <= ngroup; g++){
	      ih ++;
	      pntr_ag(f,g) = ih;
	      for(h = 1; h <= nsex; h++){
	        ig ++;
	        n_area(ig) = f;
	        n_group(ig) = g;
	        n_sex(ig) = h;
	        pntr_ags(f,g,h) = ig;
	        if(f == 1){
	          is++;
	          pntr_gs(g,h) = is;
	        }
	      }
	    }
	  }
	  if(verbose){
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| MODEL DIMENSION         |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| narea        "<<narea<<'\n';
	    LOG<<"| ngroup       "<<ngroup<<'\n';
	    LOG<<"| nsex         "<<nsex<<'\n';
	    LOG<<"| syr          "<<syr<<'\n';
	    LOG<<"| nyr          "<<nyr<<'\n';
	    LOG<<"| sage         "<<sage<<'\n';
	    LOG<<"| nage         "<<nage<<'\n';
	    LOG<<"| ngear        "<<ngear<<'\n';
	    LOG<<"| propfemale   "<<propfemale<<'\n';
	    LOG<<"| n_area       "<<n_area<<'\n';
	    LOG<<"| n_group      "<<n_group<<'\n';
	    LOG<<"| n_sex        "<<n_sex<<'\n';
	    LOG<<"| pntr_ag      "<<pntr_ag<<'\n';
	    LOG<<"| pntr_gs      "<<pntr_gs<<'\n';
	    LOG<<"| pntr_ags     "<<pntr_ags<<'\n';
	    LOG<<"| ----------------------- |\n\n";
	  }
	  // Check for dimension errors in projection control file.
	  if(pf_cntrl(1) < syr || pf_cntrl(3) < syr || pf_cntrl(5) < syr){
	    LOG<<"WARNING: start year in projection file control is"
	         " less than initial model year. Setting to syr.\n";
	    pf_cntrl(1) = syr;
	    pf_cntrl(3) = syr;
	    pf_cntrl(5) = syr;
	  }
	  if(pf_cntrl(2) > nyr || pf_cntrl(4) > nyr || pf_cntrl(6) > nyr){
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
	  dAllocation = dAllocation / sum(dAllocation);
	  for(int k = 1; k <= ngear; k++){
	    if(dAllocation(k) > 0)
	      fsh_flag(k) = 1;
	    else
	      fsh_flag(k) = 0;
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
	matrix la(1,n_ags,sage,nage); //length-at-age
	matrix wa(1,n_ags,sage,nage); //weight-at-age
	matrix ma(1,n_ags,sage,nage); //maturity-at-age
	LOC_CALCS
	  if(verbose){
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| Growth Parameters       |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| d_linf       "<<d_linf<<'\n';
	    LOG<<"| d_vonbk      "<<d_vonbk<<'\n';
	    LOG<<"| d_to         "<<d_to<<'\n';
	    LOG<<"| d_a          "<<d_a<<'\n';
	    LOG<<"| d_b          "<<d_b<<'\n';
	    LOG<<"| d_ah         "<<d_ah<<'\n';
	    LOG<<"| d_gh         "<<d_gh<<'\n';
	    LOG<<"| ----------------------- |\n\n";
	  }
	  // length & weight-at-age based on input growth pars
	  ma.initialize();
	  for(ig = 1; ig <= n_ags; ig++){
	    la(ig) = d_linf(ig) * (1. - exp(-d_vonbk(ig) * (age-d_to(ig))));
	    wa(ig) = d_a(ig) * pow(la(ig), d_b(ig));
	    h = n_sex(ig);
	    if(n_MAT == 0){
	      ma(ig) = plogis(age, d_ah(ig), d_gh(ig));
	    }else if(n_MAT > 0 && h != 2){
	      ma(ig) = d_maturityVector;
	    }
	  }
	  if(verbose){
	    //LOG<<"d_a = "<<d_a<<"\n\n";
	    //LOG<<"d_b = "<<d_b<<"\n\n";
	    //LOG<<"la = "<<la<<"\n\n";
	    //LOG<<"wa = "<<wa<<"\n\n";
	  }
	END_CALCS

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
	  if(verbose){
	    LOG<<"| Catch values ---------- |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| HEAD(dCatchData)        |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<dCatchData.sub(1,4)<<'\n';
	    LOG<<"| ----------------------- |\n\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| TAIL(dCatchData)        |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<dCatchData.sub(nCtNobs-3,nCtNobs)<<'\n';
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| All Catch               |\n";
	    LOG<<dCatchData<<"\n\n";
	  }
	  d3_Ct.initialize();
	  int k;
	  if(verbose){
	    //LOG<<"Debug reading of catch data\n";
	  }
	  for(int ii = 1; ii <= nCtNobs; ii++){
	    i = dCatchData(ii)(1); // year
	    k = dCatchData(ii)(2); // gear
	    f = dCatchData(ii)(3); // area
	    g = dCatchData(ii)(4); // group
	    h = dCatchData(ii)(5); // sex
	    if(verbose){
	      //LOG<<"Year = "<<i<<"\n";
	      //LOG<<"ii = row number = "<<ii<<"\n";
	      //LOG<<"i = year = "<<i<<"\n";
	      //LOG<<"k = gear = "<<k<<"\n";
	      //LOG<<"f = area = "<<f<<"\n";
	      //LOG<<"g = group = "<<g<<"\n";
	    }
	    if(h){
	      ig = pntr_ags(f,g,h);
	      if(verbose){
	        //LOG<<"h = sex = "<<h<<"\n";
	        //LOG<<"ig = area,group,sex = "<<ig<<"\n";
	        //LOG<<"Data row = dCatchData(ii) =  "<<dCatchData(ii)<<"\n";
	      }
	      if(ig == 1){
	        if(verbose){
	          //LOG<<"Multiplying catch by proportion female ("<<propfemale<<")\n";
	        }
	        d3_Ct(ig)(i)(k) = propfemale * dCatchData(ii)(7);
	      }else{
	        if(verbose){
	          //LOG<<"Multiplying catch by proportion male ("<<(1.0 - propfemale)<<")\n";
	        }
	        d3_Ct(ig)(i)(k) = (1.0 - propfemale) * dCatchData(ii)(7);
	      }
	      if(verbose){
	        //LOG<<"Catch value = d3_Ct(ig)(i)(k) = "<<d3_Ct(ig)(i)(k)<<"\n";
	      }
	    }else{
	      for(h = 1; h <= nsex; h++){
	        ig = pntr_ags(f,g,h);
	        if(verbose){
	          //LOG<<"h = sex = "<<h<<"\n";
	          //LOG<<"ig = area,group,sex = "<<ig<<"\n";
	          //LOG<<"Data row = dCatchData(ii) = "<<dCatchData(ii)<<"\n";
	        }
	        d3_Ct(ig)(i)(k) = 1.0 / nsex * dCatchData(ii)(7);
	        if(verbose){
	          //LOG<<"Catch value = d3_Ct("<<ig<<")("<<i<<")("<<k<<") = "<<d3_Ct(ig)(i)(k)<<"\n";
	        }
	      }
	      if(verbose){
	        //LOG<<"\n";
	      }
	    }
	  }
	  for(ig = 1; ig<= 2; ig++){
	    for(i = 1996; i <= 2021; i++){
	      //LOG<<"ig = "<<ig<<"\n";
	      //LOG<<"i (year) = "<<i<<"\n";
	      for(k = 1; k<=ngear; k++){
	        //LOG<<"k = gear = "<<k<<"\n";
	        //LOG<<"d3_Ct("<<ig<<")("<<i<<")("<<k<<") = "<<d3_Ct(ig)(i)(k)<<"\n";
	      }
	      //LOG<<"\n";
	    }
	    //LOG<<"\n";
	  }
	  //LOG<<"\n\n";
	END_CALCS
	// |---------------------------------------------------------------------------------|
	// | RELATIVE ABUNDANCE INDICIES (ragged array)
	// |---------------------------------------------------------------------------------|
	// | nItNobs        : number of independent surveys
	// | n_it_nobs      : number of survey observations
	// | n_survey_type  : 1 survey is proportional to vulnerable numbers
	// | n_survey_type  : 2 survey is proportional to vulnerable biomass
	// | n_survey_type  : 3 survey is proportional to vulnerable spawning biomass
	// | d3_survey_data : (iyr index(it) gear area group sex wt timing)
	// | it_wt          : relative weights for each relative abundance normalized to have a
	// |                  mean = 1 so rho = sig ^ 2 / (sig ^ 2 + tau ^ 2) holds true in variance pars.
	// |
	init_int nItNobs;
	ivector nSurveyIndex(1,nItNobs);
	init_ivector n_it_nobs(1,nItNobs);
	init_ivector n_survey_type(1,nItNobs);
	init_3darray d3_survey_data(1,nItNobs,1,n_it_nobs,1,8);
	matrix it_wt(1,nItNobs,1,n_it_nobs);
	matrix it_grp(1,nItNobs,1,n_it_nobs);
	LOC_CALCS
	  if(verbose){
	    LOG<<"| Survey indices -------- |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| HEAD(d3_survey_data)    |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<d3_survey_data(1).sub(1,4)<<'\n';
	    LOG<<"| ----------------------- |\n\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| TAIL(d3_survey_data)    |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<d3_survey_data(nItNobs).sub(n_it_nobs(nItNobs) - 3, n_it_nobs(nItNobs))<<'\n';
	    LOG<<"| ----------------------- |\n\n";
	  }
	  for(int k = 1; k <= nItNobs; k++){
	    it_wt(k) = column(d3_survey_data(k),7) + 1.e-30;
	    it_grp(k)= column(d3_survey_data(k),5);
	    nSurveyIndex(k) = d3_survey_data(k)(1,3);
	  }
	  for(int k = 1; k <= nItNobs; k++){
	    it_wt(k) = it_wt(k) / mean(it_wt(k));
	  }
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | AGE OR LENGTH COMPOSITION DATA (ragged object)
	// |---------------------------------------------------------------------------------|
	// | nAgears      -> number of age-composition matrixes, one for each gear.
	// | n_A_nobs     -> ivector for number of rows in age composition (A) matrix
	// | n_A_sage     -> imatrix for starting age in each row
	// | n_A_nage     -> imatrix for plus group age in each row
	// | inp_nscaler  -> effective sample size for iterative re-weighting in multinomial.
	// | icol_A       -> number of columns for each row in A.
	// | A            -> array of data (year,gear,area,group,sex|Data...)
	// | d3_A_obs     -> array of catch-age data only.
	// |
	init_int nAgears;
	init_ivector n_A_nobs(1,nAgears);
	init_ivector n_A_sage(1,nAgears);
	init_ivector n_A_nage(1,nAgears);
	init_vector inp_nscaler(1,nAgears);
	init_ivector n_ageFlag(1,nAgears);
	init_number dm_num_samp; // Number of samples for all DM inputs. Use only if dm_use_single_num_samp is set
	init_number dm_use_single_num_samp; // Use dm_num_samp for all samples, ignore individual sample data
	// The 6 in the next command is to remove the first 6 columns
	// from the age comp data because they are not the actual ages,
	// but the header data.
	init_3darray d3_A(1,nAgears,1,n_A_nobs,n_A_sage-6,n_A_nage);
	3darray d3_A_obs(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	matrix samp_sizes(1,nAgears,1,n_A_nobs)
	LOC_CALCS
	  if(n_A_nobs(nAgears) > 0){
	    LOG<<"| Age compositions ------  |\n";
	    LOG<<"| -----------------------  |\n";
	    LOG<<"| HEAD(d3_A) first 3 lines |"<<'\n';
	    LOG<<"| -----------------------  |\n";
	    LOG<<d3_A(1).sub(1,3)<<'\n';
	    LOG<<"| -----------------------  |\n";
	    LOG<<"| TAIL(d3_A) last line     |"<<'\n';
	    LOG<<"| -----------------------  |\n";
	    LOG<<d3_A(nAgears).sub(n_A_nobs(nAgears),n_A_nobs(nAgears))<<'\n';
	    LOG<<"| ----------------------- |\n";
	  }else{
	    LOG<<"| Age compositions ------ |\n";
	    LOG<<"| ----------------------- |\n";
	    LOG<<"| NO AGE OR LENGTH DATA   |\n";
	    LOG<<"| ----------------------- |\n";
	  }
	  // Set up Dirichlet Multinomial sample sizes
	  for(k = 1; k <= nAgears; k++){
	    if(dm_use_single_num_samp){
	      samp_sizes(k) = dm_num_samp;
	    }else{
	      samp_sizes(k) = column(d3_A(k), sage - 5);
	    }
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
	init_3darray d3_mean_wt_data(1,nMeanWt,1,nMeanWtNobs,1,7);

	LOC_CALCS
	  /*
	    This will determine the new dimension of d3_inp_wt_avg in case the backward
	    projection is needed required and rename nWtNobs to tmp_nWtNobs
	  */
	  projwt.initialize();
	  n_bf_wt_row.initialize();
	  tmp_nWtNobs.initialize();
	  for(int k=1; k<=nWtTab; k++){
	    tmp_nWtNobs(k) = nWtNobs(k);
	    projwt(k) = 1;
	    for(i=1; i<=nWtNobs(k); i++){
	      if(nWtNobs(k) > 0 && d3_inp_wt_avg(k)(i)(sage-5) < 0){
	        n_bf_wt_row(k)++;
	      }
	    }
	    if(n_bf_wt_row(k)>0){
	      for(int i=1; i<=n_bf_wt_row(k); i++){
	        int exp_nyr = fabs(d3_inp_wt_avg(k,i,sage-5)) - syr;
	        tmp_nWtNobs(k) += exp_nyr;
	      }
	      projwt(k) = -n_bf_wt_row(k);
	    }else if (n_bf_wt_row(k) == 0){
	      tmp_nWtNobs(k) = nWtNobs(k);
	      projwt(k) = 1;
	    }
	  }
	  sum_tmp_nWtNobs = sum(tmp_nWtNobs);
	END_CALCS

	3darray xinp_wt_avg(1,nWtTab,1,tmp_nWtNobs,sage-5,nage);
	matrix xxinp_wt_avg(1,sum_tmp_nWtNobs,sage-5,nage);

	LOC_CALCS
	  /*
	    This will redimension the d3_inp_wt_avg  according to tmp_nWtNobs and rename
	    the 3d array to xinp_wt_avg. Then the 3darray is converted to a matrix
	    xxinp_wt_avg
	  */
	  xinp_wt_avg.initialize();
	  xxinp_wt_avg.initialize();
	  for(int k=1; k<=nWtTab; k++){
	    ivector iroww(0,n_bf_wt_row(k));
	    iroww.initialize();
	    if(nWtNobs(k) > 0){
	      if(n_bf_wt_row(k) > 0){
	        for(i=1; i<=n_bf_wt_row(k); i++){
	          d3_inp_wt_avg(k,i,sage-5) = fabs(d3_inp_wt_avg(k, i, sage - 5));
	          iroww(i) = d3_inp_wt_avg(k, i, sage - 5) - syr + iroww(i - 1);
	          for(int jj=iroww(i);jj>=iroww(i-1)+1;jj--){
	            xinp_wt_avg(k)(jj)(sage-5) = syr + jj - iroww(i - 1) - 1;
	            xinp_wt_avg(k)(jj)(sage-4,nage) = d3_inp_wt_avg(k)(i)(sage-4,nage);
	          }
	        }
	        for(int jj = iroww(n_bf_wt_row(k))+1; jj <= tmp_nWtNobs(k); jj++){
	          xinp_wt_avg(k)(jj)(sage-5,nage) = d3_inp_wt_avg(k)(jj-iroww(n_bf_wt_row(k)))(sage-5,nage);
	        }
	      }else{
	        for(int jj = 1; jj <= tmp_nWtNobs(k); jj++){
	          xinp_wt_avg(k)(jj)(sage-5,nage) = d3_inp_wt_avg(k)(jj)(sage-5,nage);
	        }
	      }
	      int ttmp = sum(tmp_nWtNobs(1, k - 1));
	      int ttmp2 = sum(tmp_nWtNobs(1, k));
	      for(int jj=ttmp+1; jj<=ttmp2; jj++){
	        xxinp_wt_avg(jj)(sage-5,nage) = xinp_wt_avg(k)(jj-ttmp)(sage-5,nage);
	      }
	    }
	  }
	END_CALCS

	matrix dWt_bar(1,n_ags,sage,nage);
	3darray d3_wt_avg(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_dev(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_mat(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_len_age(1,n_ags,syr,nyr+1,sage,nage);

	LOC_CALCS
	  d3_wt_avg.initialize();
	  d3_wt_dev.initialize();
	  d3_wt_mat.initialize();
	  d3_len_age.initialize();
	  for(ig=1;ig<=n_ags;ig++){
	    for(int i = syr; i <= nyr; i++){
	      d3_wt_avg(ig)(i) = wa(ig);
	      //d3_wt_mat(ig)(i) = pow(elem_prod(ma(ig),wa(ig)),d_iscamCntrl(6));
	      d3_wt_mat(ig)(i) = elem_prod(ma(ig),wa(ig));
	      d3_len_age(ig)(i) = pow(wa(ig)/d_a(ig),1./d_b(ig));
	      // Insert calculations for ALK here.
	    }
	  }
	  int iyr;
	  for(i = 1; i <= sum_tmp_nWtNobs; i++){
	    iyr = xxinp_wt_avg(i,sage-5);
	    f = xxinp_wt_avg(i,sage-3);
	    g = xxinp_wt_avg(i,sage-2);
	    h = xxinp_wt_avg(i,sage-1);
	    if(h){
	      ig = pntr_ags(f,g,h);
	      dvector tmp = xxinp_wt_avg(i)(sage,nage);
	      ivector idx = getIndex(age,tmp);
	      for(unsigned int ii = 1; ii <= size_count(idx); ii++){
	        d3_wt_avg(ig)(iyr)(idx(ii)) = xxinp_wt_avg(i)(idx(ii));
	        d3_len_age(ig)(iyr)(idx(ii))= pow(d3_wt_avg(ig)(iyr)(idx(ii)) / d_a(ig), 1. / d_b(ig));
	      }
	      //d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
	      d3_wt_mat(ig)(iyr) = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
	    }else{
	      for(int h = 1; h <= nsex; h++){
	        ig = pntr_ags(f,g,h);
	        dvector tmp = xxinp_wt_avg(i)(sage,nage);
	        ivector idx = getIndex(age,tmp);
	        // Problem, array indexed differ, must loop over idx;
	        // d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
	        for(unsigned int ii = 1; ii <= size_count(idx); ii++){
	          d3_wt_avg(ig)(iyr)(idx(ii)) = xxinp_wt_avg(i)(idx(ii));
	          d3_len_age(ig)(iyr)(idx(ii)) = pow(d3_wt_avg(ig)(iyr)(idx(ii)) / d_a(ig), 1. / d_b(ig));
	        }
	        d3_wt_mat(ig)(iyr) = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
	      }
	    }
	  }
	  for(ig = 1; ig <= n_ags; ig++){
	    dWt_bar(ig) = colsum(d3_wt_avg(ig).sub(pf_cntrl(3), pf_cntrl(4)));
	    dWt_bar(ig) /= pf_cntrl(4) - pf_cntrl(3) + 1;
	    d3_wt_avg(ig)(nyr+1) = dWt_bar(ig);
	    d3_wt_mat(ig)(nyr+1) = elem_prod(dWt_bar(ig),ma(ig));
	    d3_len_age(ig)(nyr+1) = pow(dWt_bar(ig)/d_a(ig),1./d_b(ig));
	  }
	  // deviations in mean weight-at-age
	  for(ig = 1; ig <= n_ags; ig++){
	    dmatrix mtmp = trans( d3_wt_avg(ig) );
	    for(j = sage; j <= nage; j++){
	      //LOG<<sum(first_difference(mtmp(j)(syr,nyr)));
	      if(sum(first_difference(mtmp(j)(syr,nyr)))){
	        mtmp(j) = (mtmp(j)-mean(mtmp(j)(syr,nyr))) / sqrt(var(mtmp(j)(syr,nyr)));
	      }else{
	        mtmp(j) = 0;
	      }
	    }
	    d3_wt_dev(ig) = trans(mtmp);
	    if(min(d3_wt_avg(ig)) <= 0.000 && min(d3_wt_avg(ig)) != NA){
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
	init_int eof;
	LOC_CALCS
	  if(eof == 999){
	    LOG<<" ______________________________ \n";
	    LOG<<"|      END OF DATA SECTION     |\n";
	    LOG<<"|______________________________|\n";
	  }else{
	    LOG<<"\n *** ERROR READING DATA *** - eof value = "<<eof<<"\n\n";
	    exit(1);
	  }
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | VARIABLES FOR MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	matrix fmsy(1,ngroup,1,nfleet); //Fishing mortality rate at Fmsy
	matrix fall(1,ngroup,1,nfleet); //Fishing mortality based on dAllocation
	matrix  msy(1,ngroup,1,nfleet); //Maximum sustainable yield
	vector bmsy(1,ngroup); //Spawning biomass at MSY
	// number Umsy; //Exploitation rate at MSY
	vector age_tau2(1,nAgears); //MLE estimate of the variance for age comps

	LOC_CALCS
	  LOG<<" ___________________________________\n";
	  LOG<<"|      START OF CONTROL SECTION     |\n";
	  LOG<<"|___________________________________|\n\n";
	END_CALCS
	// |---------------------------------------------------------------------------------|
	// | CONTROL FILE
	// |---------------------------------------------------------------------------------|
	// |
	!! ad_comm::change_datafile_name(ControlFile);

	// |---------------------------------------------------------------------------------|
	// | Leading Parameters
	// |---------------------------------------------------------------------------------|
	// | npar             -> number of leading parameters
	// | ipar_vector      -> integer vector based on the number of areas groups sexes
	// | -1) log_ro       - unfished sage recruitment
	// | -2) steepness    - steepness of the stock-recruitment relationship
	// | -3) log_m_female - instantaneous natural mortality rate
	// | If nsex == 1:
	// | -4) log_avgrec   - average sage recruitment from syr+1 to nyr
	// | -5) log_recinit  - average sage recruitment for initialization
	// | -6) rho          - proportion of total variance for observation errors
	// | -7) vartheta     - total precision (1/variance)
	// | If nsex == 2:
	// | -4) log_m_male   - instantaneous natural mortality rate
	// | -5) log_avgrec   - average sage recruitment from syr+1 to nyr
	// | -6) log_recinit  - average sage recruitment for initialization
	// | -7) rho          - proportion of total variance for observation errors
	// | -8) vartheta     - total precision (1/variance)
	init_int npar;
	init_matrix theta_control(1,npar,1,7);

	vector theta_ival(1,npar);
	vector theta_lb(1,npar);
	vector theta_ub(1,npar);
	ivector theta_phz(1,npar);
	ivector theta_prior(1,npar);
	ivector ipar_vector(1,npar);
	LOC_CALCS
	  theta_ival = column(theta_control,1);
	  theta_lb = column(theta_control,2);
	  theta_ub = column(theta_control,3);
	  theta_phz = ivector(column(theta_control,4));
	  theta_prior = ivector(column(theta_control,5));
	  ipar_vector(1,2) = ngroup;
	  ipar_vector(3) = n_gs;
	  if(nsex == 2){
	    ipar_vector(4) = n_gs;
	    ipar_vector(5,6) = n_ag;
	    ipar_vector(7,8) = ngroup;
	  }else{
	    ipar_vector(4,5) = n_ag;
	    ipar_vector(6,7) = ngroup;
	  }
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Initial values          |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<theta_ival<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Lower bounds            |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<theta_lb<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Upper bounds            |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<theta_ub<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Phase                   |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<theta_phz<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Prior type              |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<theta_prior<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | CONTROLS PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears                  |
	// |---------------------------------------------------------------------------------|
	// |

	init_ivector nCompIndex(1,nAgears);
	init_ivector nCompLikelihood(1,nAgears);
	init_vector dMinP(1,nAgears);
	init_vector dEps(1,nAgears);
	init_vector nPhz_dm_vec(1,nAgears);
	init_ivector nPhz_age_tau2(1,nAgears);
	init_ivector nPhz_phi1(1,nAgears);
	init_ivector nPhz_phi2(1,nAgears);
	init_ivector nPhz_df(1,nAgears);
	number nPhz_dm;
	LOC_CALCS
	  nPhz_dm = nPhz_dm_vec(1);
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Gear indices            |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<nCompIndex<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| ----------------------- |\n";
	  LOG<<"| Likelihood types        |"<<'\n';
	  LOG<<"| ----------------------- |\n";
	  LOG<<nCompLikelihood<<'\n';
	  LOG<<"| ----------------------- |\n\n";
	  LOG<<"| -------------------------------------- |\n";
	  LOG<<"| Min props for agg and tail compression |"<<'\n';
	  LOG<<"| -------------------------------------- |\n";
	  LOG<<dMinP<<'\n';
	  LOG<<"| -------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Small consts to add to comps          |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<dEps<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phases for agecomp DM estimation      |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nPhz_dm<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phases for log_age_tau2 estimation    |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nPhz_age_tau2<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phases for phi1 estimation            |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nPhz_phi1<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phases for phi2 estimation            |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nPhz_phi2<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phases for degrees of freedom         |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nPhz_df<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  // Apply effective sample size if multinomial likelihood
	  if(n_A_nobs(nAgears) > 0){
	    for(k = 1; k <= nAgears; k++){
	      dmatrix tmp = trans(trans(d3_A(k)).sub(n_A_sage(k), n_A_nage(k)));
	      for(i = 1; i <= n_A_nobs(k); i++){
	        tmp(i) = tmp(i) / sum(tmp(i));
	        if(inp_nscaler(k) > 0 && nCompLikelihood(k) == 2){
	          tmp(i) = tmp(i) * inp_nscaler(k);
	        }
	      }
	      d3_A_obs(k) = tmp;
	    }
	  }
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | CONTROLS FOR SELECTIVITY OPTIONS
	// |---------------------------------------------------------------------------------|
	// | 12 different options for modelling selectivity which are summarized here:
	// | isel_npar  -> ivector for # of parameters for each gear.
	// | jsel_npar  -> ivector for the number of rows for time-varying selectivity.
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
	// |    13.....age-based selectivity coefficients with age_min-age_max parameters.
	// |
	// | selex_controls (1-10)
	// |  1  -> isel_type - switch for selectivity.
	// |  2  -> ahat (sel_type=1) - age-at-50% vulnerability for logistic function or (sel_type=13) -age_min
	// |  3  -> ghat (sel_type=1) - std at 50% age of vulnerability for logistic function or (sel_type=13) -age_max
	// |  4  -> age_nodes - No. of age-nodes for bicubic spline.
	// |  5  -> yr_nodes  - No. of year-nodes for bicubic spline.
	// |  6  -> sel_phz   - phase for estimating selectivity parameters.
	// |  7  -> lambda_1  - penalty weight for 2nd difference in selectivity.
	// |  8  -> lambda_2  - penalty weight for dome-shaped selectivity.
	// |  9  -> lambda_3  - penalty weight for 2nd difference in time-varying selectivity.
	// |  10 -> Number of discrete selectivity blocks.

	init_matrix selex_controls(1,12,1,ngear);
	ivector isel_npar(1,ngear);
	ivector jsel_npar(1,ngear);
	ivector isel_type(1,ngear);
	ivector sel_phz(1,ngear);
	ivector n_sel_blocks(1,ngear);
	vector ahat_agemin_f(1,ngear);
	vector ghat_agemax_f(1,ngear);
	vector ahat_agemin_m(1,ngear);
	vector ghat_agemax_m(1,ngear);
	vector age_nodes(1,ngear);
	vector yr_nodes(1,ngear);
	vector lambda_1(1,ngear);
	vector lambda_2(1,ngear);
	vector lambda_3(1,ngear);
	LOC_CALCS
	  ahat_agemin_f = selex_controls(2);
	  ghat_agemax_f = selex_controls(3);
	  ahat_agemin_m = selex_controls(4);
	  ghat_agemax_m = selex_controls(5);
	  age_nodes = selex_controls(6);
	  yr_nodes = selex_controls(7);
	  lambda_1 = selex_controls(9);
	  lambda_2 = selex_controls(10);
	  lambda_3 = selex_controls(11);
	  isel_type = ivector(selex_controls(1));
	  sel_phz = ivector(selex_controls(8));
	  n_sel_blocks = ivector(selex_controls(12));
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Selectivity types                     |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<isel_type<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Female Age/length at 50% selectivity  |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<ahat_agemin_f<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Female STD ot 50% selectivity         |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<ghat_agemax_f<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Male Age/length at 50% selectivity    |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<ahat_agemin_m<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Male STD ot 50% selectivity           |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<ghat_agemax_m<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| No. age nodes for each gear           |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<age_nodes<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| No. year nodes for 2d spline          |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<yr_nodes<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Phase of estimation                   |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<sel_phz<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Penalty wt for 2nd differences        |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<lambda_1<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Penalty wt for dome-shaped            |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<lambda_2<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Penalty wt for time-varting selex     |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<lambda_3<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Number of selectivity blocks          |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<n_sel_blocks<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	END_CALCS

	init_imatrix sel_blocks(1,ngear,1,n_sel_blocks);
	LOC_CALCS
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Selectivity blocks                    |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<sel_blocks<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  isel_npar.initialize();
	  for(i = 1; i <= ngear; i++){
	    jsel_npar(i) = 1;
	    switch(isel_type(i)){
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
	        isel_npar(i) = (ghat_agemax_f(i)-ahat_agemin_f(i)+1);
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
	// | Controls for fitting to mean weight data
	// |---------------------------------------------------------------------------------|
	// | fitMeanWt         -> 1 = fit to annual mean weights; 0 = do not fit to annual mean weights
	// | nMeanWtCV         -> Number of annual mean weight series
	// | vector weight_sig -> sd for likelihood for fitting to annual mean weight (one for each series)
	init_int fitMeanWt;
	init_int nMeanWtCV;
	init_vector weight_sig(1,nMeanWtCV);

	LOC_CALCS
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Number or surveys                     |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<nits<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Prior type for Q                      |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<q_prior<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| mu_log_q for Prior                    |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<mu_log_q<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| sd_log_q for Prior                    |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<sd_log_q<<'\n';
	  LOG<<"| ------------------------------------- |\n\n";
	END_CALCS
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
	// |      If this is greater than 0, the "slow" fmsy routine will
	// |      be used instead of the newton-rhapson routine. The MSY-
	// |      based reference points in the report file will reflect this.
	// |      If greater than 0, a file called TEST_frp.rep will be produced.
	// | 14-> number of prospective years to start estimation from syr.
	// | 15-> switch for generating selex based on IFD and cohort biomass
	// | 16-> toggle to fit to annual mean weights for commercial catch
	// | 17-> toggle to perform "slow" fmsy test (runs model out 100 years)
	// |      this produces a file called TEST_frp.rep.
	// | 18-> precision for F for the "slow" fmsy calculations (only used if
	// |      control 18 is 1). This must be a maximum of 0.0001 or the
	// |      program will stop.
	// | 19-> maximum F for the "slow" fmsy calculations (only used if
	// |      control 18 is 1). If this is greater than 1, a warning will
	// |      be issued because it will take a long time to run.
	// | 20-> Report b0 only in mcmc calculations even if the "slow" msy
	// |      routine was run. MSY-based reference points will not be
	// |      output for MCMCs. This control only matters if control 13
	// |      is greater than 0.

	init_vector d_iscamCntrl(1,21);
	int verbose;
	init_int eofc;
	LOC_CALCS
	  LOG<<"| ------------------------------------- |\n";
	  LOG<<"| Micellaneous controls                 |"<<'\n';
	  LOG<<"| ------------------------------------- |\n";
	  for(unsigned int i = 1; i <= d_iscamCntrl.size(); i++){
	    LOG<<"#"<<i<<" = "<<d_iscamCntrl(i)<<'\n';
	  }
	  LOG<<"| ------------------------------------- |\n";
	  for(int ig = 1; ig <= n_ags; ig++){
	    for(int i = syr; i <= nyr; i++){
	      d3_wt_mat(ig)(i) = pow(d3_wt_mat(ig)(i),d_iscamCntrl(6));
	    }
	  }
	  if(eofc == 999){
	    LOG<<" ______________________________ \n";
	    LOG<<"|     END OF CONTROL FILE      |\n";
	    LOG<<"|______________________________|\n";
	  }else{
	    LOG<<"\n ***** ERROR READING CONTROL FILE ***** \n";
	    exit(1);
	  }
	  if((d_iscamCntrl(13) || d_iscamCntrl(18)) && d_iscamCntrl(18) > 0.0001){
	    cerr<<"Error - you have set the precision for the slow msy calculations"
	    " too high. The maximum is 0.0001.\n";
	    exit(1);
	  }
	  if((d_iscamCntrl(13) || d_iscamCntrl(18)) && d_iscamCntrl(18) < 0.000001){
	    cout<<"Warning - you have set the precision for the slow msy"
	    " below 0.000001. This may cause the program to run for a longer time.\n";
	  }
	  if((d_iscamCntrl(13) || d_iscamCntrl(18)) && d_iscamCntrl(19) > 1){
	    cout<<"Warning - you have set the maximum F value above 1. This may cause"
	    " the program to run for a longer time.\n";
	  }
	  verbose = d_iscamCntrl(1);
	END_CALCS

	int nf;

	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENSIONS FOR NEGATIVE LOG LIKELIHOODS
	// |---------------------------------------------------------------------------------|
	// | ilvec[1,6,7,8]      -> number of fishing gears (ngear)
	// | ilvec[2]            -> number of surveys       (nItNobs)
	// | ilvec[3]            -> number of age-compisition data sets (nAgears)
	// | ilvec[4]            -> container for recruitment deviations.
	// | ilvec[5]            -> number of annual mean weight datasets.
	ivector ilvec(1,8);
	!! ilvec = ngear;
	!! ilvec(1) = 1;
	!! ilvec(2) = nItNobs;
	!! ilvec(3) = nAgears;
	!! ilvec(4) = ngroup;
	!! ilvec(5) = nMeanWt;

	// |---------------------------------------------------------------------------------|
	// | RETROSPECTIVE ADJUSTMENT TO nyrs
	// |---------------------------------------------------------------------------------|
	// | - Don\'t read any more input data from here on in.
	// | - Modifying nyr to allow for retrospective analysis.
	// | - If retro_yrs > 0, then ensure that pf_cntrl arrays are not greater than nyr,
	// |   otherwise arrays for mbar will go out of bounds.
	// | - Reduce ft_count so as not to bias estimates of ft.
	// | - Establish retrospective counter for Composition data n_naa;

	ivector n_naa(1,nAgears);
	!! nyr = nyr - retro_yrs;
	LOC_CALCS
	  if(retro_yrs){
	    if(pf_cntrl(2) > nyr)
	      pf_cntrl(2) = nyr;
	    if(pf_cntrl(4) > nyr)
	      pf_cntrl(4) = nyr;
	    if(pf_cntrl(6) > nyr)
	      pf_cntrl(6) = nyr;
	  }
	  for( i = 1; i <= nCtNobs; i++ ){
	    if( dCatchData(i)(1) > nyr )
	      ft_count --;
	  }
	  // Retrospective counter for n_A_nobs
	  n_naa.initialize();
	  for(k = 1; k <= nAgears; k++){
	    for(i = 1; i <= n_A_nobs(k); i++){
	      iyr = d3_A(k)(i)(n_A_sage(k) - 6);	//index for year
	      if(iyr <= nyr)
	        n_naa(k)++;
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
	  if(pf_cntrl(1) < syr)
	    pf_cntrl(1) = syr;
	  if(pf_cntrl(3) < syr)
	    pf_cntrl(3) = syr;
	  if(pf_cntrl(5) < syr)
	    pf_cntrl(5) = syr;
	  for(i = 1; i <= nCtNobs; i++){
	    if(dCatchData(i)(1) < syr )
	      ft_count --;
	  }
	  // Prospective counter for n_A_nobs
	  n_saa.initialize();
	  n_saa = 1;
	  for(int k = 1; k <= nAgears; k++){
	    for(int i = 1; i <= n_A_nobs(k); i++){
	      iyr = d3_A(k)(i)(n_A_sage(k) - 6); //index for year
	      if(iyr < syr){
	        n_saa(k)++;
	      }
	    }
	    }
	END_CALCS

	friend_class OperatingModel;

INITIALIZATION_SECTION
	theta theta_ival;
	phi1 0.01;

PARAMETER_SECTION
	// |---------------------------------------------------------------------------------|
	// | LEADING PARAMTERS
	// |---------------------------------------------------------------------------------|
	// | Initialized in the INITIALIZATION_SECTION with theta_ival from control file.
	// | Change to init_bounded_vector_vector.
	// | theta[1] -> log_ro, or log_msy
	// | theta[2] -> steepness(h), or log_fmsy
	// | theta[3] -> log_m_female
	// | if nsex == 1:
	// | theta[4] -> log_avgrec
	// | theta[5] -> log_recinit
	// | theta[6] -> rho
	// | theta[7] -> vartheta
	// | if nsex == 2:
	// | theta[4] -> log_m_male
	// | theta[5] -> log_avgrec
	// | theta[6] -> log_recinit
	// | theta[7] -> rho
	// | theta[8] -> vartheta
	init_bounded_vector_vector theta(1,npar,1,ipar_vector,theta_lb,theta_ub,theta_phz);

	// |---------------------------------------------------------------------------------|
	// | DIRICHLET MULTINOMIAL PARAMETERS
	// |---------------------------------------------------------------------------------|
	// Set up Dirichlet Multinomial parameters
	init_bounded_vector log_phi(1,nAgears,0,10,nPhz_dm);
	vector temp_n(1,nage)

	// |---------------------------------------------------------------------------------|
	// | SELECTIVITY PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | This is a bounded matrix vector where the dimensions are defined by the
	// |   selectivity options specified in the control file.
	// | There are 1:ngear arrays, having jsel_npar rows and isel_npar columns.
	// | If the user has not specified -ainp or -binp, the initial values are set
	// |   based on ahat and ghat in the control file for logistic selectivities.
	// | NB sel_par is in log space.
	init_bounded_matrix_vector sel_par_f(1,ngear,1,jsel_npar,1,isel_npar,-25.,25.,sel_phz);
	init_bounded_matrix_vector sel_par_m(1,ngear,1,jsel_npar,1,isel_npar,-25.,25.,sel_phz);

	LOC_CALCS
	  //LOG<<"global_parfile = "<<global_parfile<<"\n\n";
	  if(!global_parfile){
	    for(int k = 1; k <= ngear; k++){
	      if(isel_type(k) == 1 ||
	         isel_type(k) == 6 ||
	        (isel_type(k) >= 7 &&
	         isel_type(k) <= 12)){
	        if(verbose){
	          //LOG<<"SEL_PARS: gear = "<<k<<", n_sel_blocks(k) = "<<n_sel_blocks(k)<<"\n";
	        }
	        for(int j = 1; j <= n_sel_blocks(k); j++ ){
	          double uu = 0;
	          sel_par_f(k,j,1) = log(ahat_agemin_f(k) * exp(uu));
	          sel_par_f(k,j,2) = log(ghat_agemax_f(k));
	          if(nsex == 2){
	            sel_par_m(k,j,1) = log(ahat_agemin_m(k) * exp(uu));
	            sel_par_m(k,j,2) = log(ghat_agemax_m(k));
	          }
	        }
	      }else if(isel_type(k) == 13){
	        for(int j = 1; j <= n_sel_blocks(k); j++){
	          double dd = 1.e-8;
	          double stp = 1.0 / (ghat_agemax_f(k) - ahat_agemin_f(k));
	          sel_par_f(k)(j).fill_seqadd(dd,stp);
	          if(nsex == 2){
	            stp = 1.0 / (ghat_agemax_m(k) - ahat_agemin_m(k));
	            sel_par_m(k)(j).fill_seqadd(dd,stp);
	          }
	        }
	      }
	    }
	  }
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY RATE PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | Estimate all fishing mortality rates in log-space.

	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,1);
	LOC_CALCS
	  log_ft_pars = log(0.10);
	  //LOG<<"ft_count = "<<ft_count<<"\n";
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | INITIAL AND ANNUAL RECRUITMENT
	// |---------------------------------------------------------------------------------|
	// | Estimate single mean initial recruitment and deviations for each initial
	// |   cohort from sage+1 to nage. (Rinit + init_log_rec_devs)
	// | Estimate mean overal recruitment and annual deviations from syr to nyr.
	// | d_iscamCntrl(5) is a flag to initialize the model at unfished recruitment (ro),
	// |   if this is true, then do not estimate init_log_rec_devs
	// | [ ] - TODO add dev contstraint for rec_devs in calc_objective_function.

	!! int init_dev_phz = 2;
	!! if(d_iscamCntrl(5)) init_dev_phz = -1;
	// init_log_rec_devs is small omega in the model documentation
	init_bounded_matrix init_log_rec_devs(1,n_ag,sage+1,nage,-15.,15.,init_dev_phz);

	// log_rec_devs is small phi in the model documentation
	init_bounded_matrix log_rec_devs(1,n_ag,syr,nyr,-15.,15.,2);

	// |---------------------------------------------------------------------------------|
	// | DEVIATIONS FOR NATURAL MORTALITY BASED ON CUBIC SPLINE INTERPOLATION
	// |---------------------------------------------------------------------------------|
	// | Estimating trends in natural mortality rates, where the user specified the
	// |   number of knots (d_iscamCntrl(12)) and the std in M in the control file, and the phase
	// |   in which to estimate natural mortality devs (d_iscamCntrl(10)).  If the phase is neg.
	// |   then natural mortality rate deviates are not estimated and M is assumed const.
	// | This model is implemented as a random walk, where M{t+1} = M{t} + dev.

	!! int m_dev_phz = -1;
	!! m_dev_phz = d_iscamCntrl(10);
	!! int n_m_devs = d_iscamCntrl(12);
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
	// | the value that ADMB will minimize, called objfun in iSCAM
	objective_function_value objfun;

	// |---------------------------------------------------------------------------------|
	// | POPULATION VARIABLES
	// |---------------------------------------------------------------------------------|
	// | m_bar : Average natural mortality rate from syr to nyr.
	number m_bar;

	// |---------------------------------------------------------------------------------|
	// | POPULATION VECTORS
	// |---------------------------------------------------------------------------------|
	// | ro          -> theoretical unfished age-sage recruits.
	// | bo          -> theoretical unfished spawning biomass (MSY-based ref point).
	// | sbo         -> unfished spawning biomass at the time of spawning.
	// | kappa       -> Goodyear recruitment compensation ratio K = 4h/(1-h); h=K/(4+K)
	// | so          -> Initial slope (max R/S) of the stock-recruitment relationship.
	// | beta        -> Density dependent term in the stock-recruitment relationship.
	// | m           -> Instantaneous natural mortality rate by nsex
	// | log_avgrec  -> Average sage recruitment(syr-nyr,area,group).
	// | log_recinit -> Avg. initial recruitment for initial year cohorts(area,group).
	// | log_m_devs  -> annual deviations in natural mortality.
	// | q           -> conditional MLE estimates of q in It=q*Bt*exp(epsilon)
	// | ct          -> predicted catch for each catch observation
	// | eta         -> standardized log residual (log(obs_ct)-log(ct))/sigma_{ct}
	// | rho         -> Proportion of total variance associated with obs error.
	// | varphi      -> Total precision of CPUE and Recruitment deviations.
	// | sig         -> STD of the observation errors in relative abundance data.
	// | tau         -> STD of the process errors (recruitment deviations).

	vector ro(1,ngroup);
	vector bo(1,ngroup);
	vector no(1,ngroup);
	vector sbo(1,ngroup);
	vector kappa(1,ngroup);
	vector steepness(1,ngroup);
	vector so(1,ngroup);
	vector beta(1,ngroup);
	vector m(1,n_gs);
	vector log_avgrec(1,n_ag);
	vector log_recinit(1,n_ag);
	vector q(1,nItNobs);
	vector ct(1,nCtNobs);
	vector eta(1,nCtNobs);
	vector log_m_devs(syr+1,nyr);
	vector rho(1,ngroup);
	vector varphi(1,ngroup);
	vector sig(1,ngroup);
	vector tau(1,ngroup);

	// |---------------------------------------------------------------------------------|
	// | MATRIX OBJECTS
	// |---------------------------------------------------------------------------------|
	// | log_rt              -> age-sage recruitment for initial years and annual recruitment.
	// | catch_df            -> Catch data_frame (year,gear,area,group,sex,type,obs,pred,resid)
	// | eta                 -> log residuals between observed and predicted total catch.
	// | nlvec               -> matrix for negative loglikelihoods.
	// | epsilon             -> residuals for survey abundance index
	// | it_hat              -> predicted survey index (no need to be differentiable)
	// | qt                  -> catchability coefficients (time-varying)
	// | sbt                 -> spawning stock biomass by group used in S-R relationship.
	// | bt                  -> average biomass by group used for stock projection
	// | rt                  -> predicted sage-recruits based on S-R relationship.
	// | delta               -> residuals between estimated R and R from S-R curve (process err)
	// | annual_mean_weight  -> ragged matrix of estimated annual mean weights for each gear with empirical annual mean weight observations

	matrix log_rt(1,n_ag,syr-nage+sage,nyr);
	matrix nlvec(1,8,1,ilvec);
	matrix epsilon(1,nItNobs,1,n_it_nobs);
	matrix std_epsilon(1,nItNobs,1,n_it_nobs);
	matrix it_hat(1,nItNobs,1,n_it_nobs);
	matrix qt(1,nItNobs,1,n_it_nobs);
	matrix sbt(1,ngroup,syr,nyr+1);
	matrix bt(1,ngroup,syr,nyr+1);
	matrix rt(1,ngroup,syr+sage,nyr);
	matrix delta(1,ngroup,syr+sage,nyr);
	matrix annual_mean_weight(1,nMeanWt,1,nMeanWtNobs);
	matrix obs_annual_mean_weight(1,nMeanWt,1,nMeanWtNobs);

	// |---------------------------------------------------------------------------------|
	// | THREE DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | ft         -> Mean fishing mortality rates for (area-sex, gear, year)
	// | F          -> Instantaneous fishing mortality rate for (group,year,age)
	// | M          -> Instantaneous natural mortality rate for (group,year,age)
	// | Z          -> Instantaneous total mortality rate Z = M + F for (group,year,age)
	// | S          -> Annual survival rate exp(-Z) for (group,year,age)
	// | N          -> Numbers-at-age for (group,year+1,age)
	// | A_hat      -> ragged matrix for predicted age-composition data.
	// | A_nu       -> ragged matrix for age-composition residuals.
	// | vbt        -> vulnerable biomass to all gears
	// |
	3darray ft(1,n_ags,1,ngear,syr,nyr);
	3darray F(1,n_ags,syr,nyr,sage,nage);
	3darray M(1,n_ags,syr,nyr,sage,nage);
	3darray Z(1,n_ags,syr,nyr,sage,nage);
	3darray S(1,n_ags,syr,nyr,sage,nage);
	3darray N(1,n_ags,syr,nyr+1,sage,nage);
	3darray A_hat(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	3darray A_nu(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	3darray vbt(1,ngroup,1,ngear,syr,nyr+1);

	// |---------------------------------------------------------------------------------|
	// | FOUR DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | log_sel    -> Selectivity for (gear, group, year, age)
	// | Chat       -> Predicted catch-age array for (gear, group, year, age)
	// |

	4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);

	// |---------------------------------------------------------------------------------|
	// | SDREPORT VARIABLES AND VECTORS
	// |---------------------------------------------------------------------------------|
	// | sd_depletion -> Predicted spawning biomass depletion level bt/Bo
	// | sd_log_sbt   -> Log Spawning biomass for each group.
	// |
	sdreport_vector sd_depletion(1,ngroup);
	sdreport_matrix sd_log_sbt(1,ngroup,syr,nyr+1);

	vector snat(1,n_gs); //natural survival
	vector sfished(1,n_ags); //natural survival
	vector wbar(1,n_gs);
	matrix annual_mean_wt(1,n_ags,syr,nyr);

PRELIMINARY_CALCS_SECTION
	// |---------------------------------------------------------------------------------|
	// | Run the model with input parameters to simulate real data.
	// |---------------------------------------------------------------------------------|
	// | - nf is a function evaluation counter.
	//LOG<<"PRELIMINARY_CALCS_SECTION\n";
	nf = 0;
	if(testMSY){
	  testMSYxls();
	}

RUNTIME_SECTION
  maximum_function_evaluations 2000, 2000, 2000, 25000, 25000
  convergence_criteria 0.01, 0.01, 1.e-3, 1.e-4, 1.e-5

PROCEDURE_SECTION
	//LOG<<"PROCEDURE_SECTION\n";
	if(d_iscamCntrl(5) == 2)
	  d_iscamCntrl(5) = 0; //This control determines whether population is unfished in syr (0=false)
	initParameters();
	calcSelectivities(isel_type);
	calcTotalMortality();
	calcNumbersAtAge();
	//LOG<<"After NumbersAtAge()\n";
	calcTotalCatch();
	//LOG<<"After calcTotalCatch()\n";
	calcComposition();
	//LOG<<"After calcComposition()\n";
	calcSurveyObservations();
	//LOG<<"After calcSurveyObservations()\n";
	calcStockRecruitment();
	//LOG<<"After calcStockRecruitment()\n";
	calcAnnualMeanWeight();
	//LOG<<"After calcAnnualMeanWeight()\n";
	calcObjectiveFunction();
	//LOG<<"After calcObjectiveFunction()\n";
	if(sd_phase()){
	  calcSdreportVariables();
	}
	if(mc_phase()){
	  mcmcPhase = 1;
	}
	if(mceval_phase()){
	  LOG<<"Running mceval phase\n";
	  mcmc_output();
	}

FUNCTION void calcSdreportVariables()
	sd_depletion.initialize();
	sd_log_sbt.initialize();
	for(g = 1 ; g <= ngroup; g++){
	  sd_depletion(g) = sbt(g)(nyr) / sbo(g);
	  sd_log_sbt(g) = log(sbt(g));
	}

	/*
	  Purpose: This function extracts the specific parameter values from the theta vector
	    to initialize the leading parameters in the model.
	  NOTES:
	    Variance partitioning:
	    Estimating total variance as = 1 / precision
	    and partition variance by rho = sig ^ 2 / (sig ^ 2 + tau ^ 2).
	    E.g. if sig = 0.2 and tau = 1.12 then
	    rho = 0.2 ^ 2 / (0.2 ^ 2 + 1.12 ^ 2) = 0.03090235
	    the total variance is kappa^2 = sig ^ 2 + tau ^ 2 = 1.2944
	  TODO list:
	    [ ] - Alternative parameterization using MSY and FMSY as leading parameters (Martell).
	    [*] - avg recruitment limited to area, may consider ragged object for area & stock.
	*/
FUNCTION void initParameters()
	int ih;
	ro = mfexp(theta(1));
	steepness = theta(2);
	m(1) = exp(theta(3,1));
	if(nsex == 2){
	  m(2) = exp(theta(4,1));
	}
	rho = theta(6 + nsex - 1);
	varphi = sqrt(1.0 / theta(7 + nsex - 1));
	sig = elem_prod(sqrt(rho) , varphi);
	tau = elem_prod(sqrt(1.0 - rho) , varphi);
	for(ih = 1; ih <= n_ag; ih++){
	  log_avgrec(ih)  = theta(4 + nsex - 1, ih);
	  log_recinit(ih) = theta(5 + nsex - 1, ih);
	}
	switch(int(d_iscamCntrl(2))){
	  case 1:
	    // Beverton-Holt model
	    kappa = elem_div(4. * steepness, (1. - steepness));
	    break;
	  case 2:
	    // Ricker model
	    kappa = pow((5. * steepness), 1.25);
	    break;
	}

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
	RETURN_ARRAYS_INCREMENT();
	int nodes = size_count(spline_coffs);
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	//some testing here
	/*
	  dvar_vector spline_nodes(1,nodes);
	  spline_nodes.fill_seqadd(-0.5,1./(nodes-1));
	  LOG<<spline_nodes<<'\n';
	  vcubic_spline_function test_ffa(ia,spline_nodes);
	  LOG<<test_ffa(fa)<<'\n';
	  exit(1);
	*/
	return(ffa(fa));

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs, const dvector& la)
	// Interplolation for length-based selectivity coefficeients
	RETURN_ARRAYS_INCREMENT();
	int nodes = size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(ffa(fa));

	/*
	  * @brief cubic spline interpolation
	  * @details Uses cubic spline interpolatoin for data type variables based on a
	  * vector of spline coefficients, or nodes, and independent points.
	  * The nodes are rescaled to 0-1.  This function does not extrapolate beyond the
	  * independent points.
	  * @param spline_coffs a data vector of spline coefficients (nodes)
	  * @param la a vector of independent points for use in interpolation.
	  * @return A data vector containing the interpolated points.
	*/
FUNCTION dvector cubic_spline(const dvector& spline_coffs, const dvector& la)
	// Interplolation for length-based selectivity coefficeients
	//RETURN_ARRAYS_INCREMENT();
	int nodes = size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	//RETURN_ARRAYS_DECREMENT();
	return(value(ffa(fa)));

	/*
	  Purpose: This function loops over each of ngears and calculates the corresponding
	  selectivity coefficients for that gear in each year.  It uses a switch
	  statement based on isel_type to determine which selectivty function to use
	  for each particular gear that is specified in the control file.  See NOTES
	  below for more information on selectivity models.
	  Arguments:
	    isel_type -> an ivector with integers that determine what selectivity model to use.

	  NOTES:
	    The following is a list of the current selectivity models that are implemented:
	    1) Logistic selectivity with 2 parameters.
	    2) Age-specific selectivity coefficients with (nage-sage) parameters.
	       and the last two age-classes are assumed to have the same selectivity.
	    3) A reduced age-specific parameter set based on a bicubic spline.
	    4) Time varying cubic spline.
	    5) Time varying bicubic spline (2d version).
	    6) Fixed logistic.
	    7) Logistic selectivity based on relative changes in mean weight at age
	    8) Time varying selectivity based on logistic with deviations in
	       weights at age (3 estimated parameters).
	   11) Logistic selectivity with 2 parameters based on mean length.
	   12) Length-based selectivity using cubic spline interpolation.
	   13) The bicubic_spline function is located in statsLib.h
	  TODO list:
	    [*] add an option for length-based selectivity.  Use inverse of
	        allometric relationship w_a = a*l_a^b; to get mean length-at-age from
	        empirical weight-at-age data, then calculate selectivity based on
	        mean length. IMPLEMENTED IN CASE 11
	    [*] change index for gear loop from j to k, and be consistent with year (i) and
	        age (j), and sex (h) indexing.
	*/
FUNCTION void calcSelectivities(const ivector& isel_type)
	int ig, i, j, k, byr, bpar, kgear;
	double tiny = 1.e-10;
	dvariable p1, p2, p3;
	dvar_vector age_dev = age;
	dvar_matrix t1;
	dvar_matrix tmp(syr, nyr - 1, sage, nage);
	dvar_matrix tmp2(syr, nyr, sage, nage);
	dvar_matrix ttmp2(sage, nage, syr, nyr);
	// Selex cSelex(age);
	// logistic_selectivity cLogisticSelex(age);
	log_sel.initialize();
	for(kgear = 1; kgear <= ngear; kgear++){
	  // The following is used to mirror another gear-type
	  // based on the absolute value of sel_phz.
	  k = kgear;
	  if(sel_phz(k) < 0){
	    k = abs(sel_phz(kgear));
	    sel_par_f(kgear) = sel_par_f(k);
	    if(nsex == 2){
	      sel_par_m(kgear) = sel_par_m(k);
	    }
	  }
	  for(ig = 1; ig <= n_ags; ig++){
	    tmp.initialize();
	    tmp2.initialize();
	    dvector iy(1,yr_nodes(k));
	    dvector ia(1,age_nodes(k));
	    byr = 1;
	    bpar = 0;
	    switch(isel_type(k)){
	      case 1: //logistic selectivity (2 parameters)
	        for(i = syr; i <= nyr; i++){
	          if(i == sel_blocks(k, byr)){
	            bpar++;
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          // LOG<<"Testing selex class"<<'\n';
	          // log_sel(k)(ig)(i) = log(cSelex.logistic(sel_par(k)(bpar)));
	          // log_sel(k)(ig)(i) = log(cLogisticSelex(sel_par(k)(bpar)));
	          // log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny );
	          if(ig == 1){
	            p1 = mfexp(sel_par_f(k, bpar, 1));
	            p2 = mfexp(sel_par_f(k, bpar, 2));
	          }else{
	            p1 = mfexp(sel_par_m(k, bpar, 1));
	            p2 = mfexp(sel_par_m(k, bpar, 2));
	          }
	          log_sel(kgear)(ig)(i) = log(plogis(age, p1, p2) + tiny);
	          if(verbose & (i == nyr)){
	            //LOG<<"sel_par_f(k, bpar, 1) = "<<sel_par_f(k, bpar, 1)<<"\n";
	            //LOG<<"sel_par_f(k, bpar, 2) = "<<sel_par_f(k, bpar, 2)<<"\n";
	            //LOG<<"sel_par_m(k, bpar, 1) = "<<sel_par_m(k, bpar, 1)<<"\n";
	            //LOG<<"sel_par_m(k, bpar, 2) = "<<sel_par_m(k, bpar, 2)<<"\n";
	            //LOG<<"k = "<<k<<", byr = "<<byr<<"\nsel_blocks(k, byr)\n"<<sel_blocks(k, byr)<<"\n";
	            //LOG<<"kgear = "<<kgear<<", ig == "<<ig<<", i = "<<i<<"\n";
	            //LOG<<"age = "<<age<<", p1 = "<<p1<<", p2 = "<<p2<<"\n";
	            //LOG<<"log_sel(kgear)(ig)(i) = "<<log_sel(kgear)(ig)(i)<<"\n\n";
	          }
	        }
	        break;
	      case 6: // fixed logistic selectivity
	        if(ig == 1){
	          p1 = mfexp(sel_par_f(k, 1, 1));
	          p2 = mfexp(sel_par_f(k, 1, 2));
	        }else{
	          p1 = mfexp(sel_par_m(k, 1, 1));
	          p2 = mfexp(sel_par_m(k, 1, 2));
		}
	        for(i = syr; i <= nyr; i++){
	          log_sel(kgear)(ig)(i) = log(plogis(age, p1, p2));
	        }
	        break;
	      case 2: // age-specific selectivity coefficients
	        for(i=syr; i<=nyr; i++){
	          if(i == sel_blocks(k,byr)){
	            bpar ++;
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          for(j = sage; j <= nage - 1; j++){
	            log_sel(kgear)(ig)(i)(j) = sel_par_f(k)(bpar)(j-sage+1);
	          }
	          log_sel(kgear)(ig)(i,nage) = log_sel(kgear)(ig)(i,nage-1);
	        }
	        break;
	      case 3: // cubic spline
	        for(i = syr; i < nyr; i++){
	          if(i == sel_blocks(k,byr)){
	            bpar ++;
	            log_sel(k)(ig)(i) = cubic_spline(sel_par_f(k)(bpar));
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          log_sel(kgear)(ig)(i+1) = log_sel(k)(ig)(i);
	        }
	        break;
	      case 4: // time-varying cubic spline every year
	        for(i = syr; i <= nyr; i++){
		  if(ig == 1){
	            log_sel(kgear)(ig)(i) = cubic_spline(sel_par_f(k)(i-syr+1));
		  }else if(ig == 2){
	            log_sel(kgear)(ig)(i) = cubic_spline(sel_par_m(k)(i-syr+1));
		  }
	        }
	        break;
	      case 5: // time-varying bicubic spline
	        ia.fill_seqadd(0,1./(age_nodes(k)-1));
	        iy.fill_seqadd( 0,1. / (yr_nodes(k)-1));
	        bicubic_spline(iy,ia,sel_par_f(k),tmp2);
	        log_sel(kgear)(ig) = tmp2;
	        break;
	      case 7:
	        // time-varying selectivity based on deviations in weight-at-age
	        // CHANGED This is not working and should not be used. (May 5, 2011)
	        // SkDM:  I was not able to get this to run very well.
	        // AUG 5, CHANGED so it no longer has the random walk component.
	        p1 = mfexp(sel_par_f(k,1,1));
	        p2 = mfexp(sel_par_f(k,1,2));
	        for(i = syr; i <= nyr; i++){
	          dvar_vector tmpwt = log(d3_wt_avg(ig)(i) * 1000) / mean(log(d3_wt_avg(ig) * 1000.));
	          log_sel(kgear)(ig)(i) = log( plogis(tmpwt,p1,p2)+tiny );
	        }
	        break;
	      case 8:
	        //Alternative time-varying selectivity based on weight
	        //deviations (d3_wt_dev) d3_wt_dev is a matrix(syr,nyr+1,sage,nage)
	        //p3 is the coefficient that describes variation in log_sel.
	        p1 = mfexp(sel_par_f(k, 1, 1));
	        p2 = mfexp(sel_par_f(k, 1, 2));
	        p3 = sel_par_f(k, 1, 3);
	        for(i = syr; i <= nyr; i++){
	          tmp2(i) = p3 * d3_wt_dev(ig)(i);
	          //log_sel(kgear)(ig)(i) = log(plogis<dvar_vector>(age, p1, p2) + tiny ) + tmp2(i);
	          log_sel(kgear)(ig)(i) = log(plogis(age, p1, p2) + tiny) + tmp2(i);
	        }
	        break;
	      case 11: // logistic selectivity based on mean length-at-age
	        for(i = syr; i <= nyr; i++){
	          if(i == sel_blocks(k,byr)){
	            bpar ++;
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          p1 = mfexp(sel_par_f(k, bpar, 1));
	          p2 = mfexp(sel_par_f(k, bpar, 2));
	          dvector len = pow(d3_wt_avg(ig)(i) / d_a(ig), 1. / d_b(ig));
	          log_sel(kgear)(ig)(i) = log(plogis(len, p1, p2));
	        }
	        break;
	      case 12: // cubic spline length-based coefficients.
	        for(i = syr; i <= nyr; i++){
	          if(i == sel_blocks(k,byr)){
	            bpar ++;
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          dvector len = pow(d3_wt_avg(ig)(i) / d_a(ig), 1. / d_b(ig));
	          log_sel(kgear)(ig)(i) = cubic_spline(sel_par_f(k)(bpar), len);
	        }
	        break;
	      case 13: // truncated age-specific selectivity coefficients
	        for(i = syr; i <= nyr; i++){
	          if(i == sel_blocks(k,byr)){
	            bpar ++;
	            if(byr < n_sel_blocks(k))
	              byr++;
	          }
	          for(j = ahat_agemin_f(k); j <= ghat_agemax_f(k); j++){
	            log_sel(k)(ig)(i)(j) = sel_par_f(k)(bpar)(j-ahat_agemin_f(k)+1);
	          }
	          for(j = ghat_agemax_f(k) + 1; j <= nage; j++){
	            log_sel(kgear)(ig)(i,j) = log_sel(kgear)(ig)(i,ghat_agemax_f(k));
	          }
	          for(j = sage; j < ahat_agemin_f(k); j++){
	              log_sel(kgear)(ig)(i,j) = log_sel(kgear)(ig)(i,ahat_agemin_f(k));
	          }
	        }
	        break;
	      default:
	        log_sel(kgear)(ig) = 0;
	        break;
	    }
	    // subtract mean to ensure mean(exp(log_sel)) == 1
	  }
	}

	/*
	  Purpose: This function calculates fishing mortality, total mortality and annual
	    surivival rates S=exp(-Z) for each age and year based on fishing mortality
	    and selectivity coefficients. Z also is updated with time-varying
	    natural mortality rates if specificed by user.
	  NOTES:
	-   Jan 5, 2012 Added catch_type to allow for catch in numbers, weight or spawn.
	    In the case of spawn on kelp (roe fisheries), the Fishing mortality does not
	      occur on the adult component.
	    Added if(catch_type(k)!=3) //exclude roe fisheries
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
	int ig, ii, i, k, l;
	int ft_counter = 0;
	dvariable ftmp;
	F.initialize();
	ft.initialize();
	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	for(ig = 1; ig <= nCtNobs; ig++){
	  i = dCatchData(ig)(1); //year
	  k = dCatchData(ig)(2); //gear
	  f = dCatchData(ig)(3); //area
	  g = dCatchData(ig)(4); //group
	  h = dCatchData(ig)(5); //sex
	  l = dCatchData(ig)(6); //type
	  if(i < syr || i > nyr)
	    continue;
	  ft_counter++;
	  if(h){
	    ii = pntr_ags(f, g, h);
	    ftmp = mfexp(log_ft_pars(ft_counter));
	    ft(ii)(k, i) = ftmp;
	    if(l != 3){
	      F(ii)(i) += ftmp * mfexp(log_sel(k)(ii)(i));
	    }
	  }else{
	    for(h = 1; h <= nsex; h++){
	      ii = pntr_ags(f, g, h);
	      ftmp = mfexp(log_ft_pars(ft_counter));
	      ft(ii)(k, i) = ftmp;
	      if(l != 3){
	        F(ii)(i) += ftmp * mfexp(log_sel(k)(ii)(i));
	      }
	      //if(last_phase() & (i == 2020)){
	      if(std::isnan(value(log_ft_pars(ft_counter)))){
	        //LOG<<"In CalcTotalMortality()\n";
	        //LOG<<"dCatchData("<<ig<<") = "<<dCatchData(ig)<<"\n";
	        //LOG<<"log_ft_pars("<<ft_counter<<") = "<<log_ft_pars(ft_counter)<<"\n";
	        //LOG<<"ftmp = "<<ftmp<<"\n";
	        //LOG<<"log_ft_pars = "<<log_ft_pars<<"\n";
	        //LOG<<"ft("<<ii<<")("<<k<<", "<<i<<") = "<<ft(ii)(k, i)<<"\n";
	        //LOG<<"F("<<ii<<")("<<i<<") = "<<F(ii)(i)<<"\n";
	        //LOG<<"mfexp(log_sel("<<k<<")("<<ii<<")("<<i<<")) = "<<mfexp(log_sel(k)(ii)(i))<<"\n";
	        //LOG<<"ft_counter = "<<ft_counter<<", l = "<<l<<"\n\n";
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
	for(ig = 1; ig <= n_ags; ig++){
	  g = n_group(ig);
	  h = n_sex(ig);
	  M(ig) = m(pntr_gs(g, h));
	  if(active(log_m_nodes)){
	    int nodes = size_count(log_m_nodes);
	    dvector im(1, nodes);
	    dvector fm(syr + 1, nyr);
	    im.fill_seqadd(0, 1. / (nodes - 1));
	    fm.fill_seqadd(0, 1. / (nyr - syr));
	    vcubic_spline_function m_spline(im, log_m_nodes);
	    log_m_devs = m_spline(fm);
	  }
	  for(i = syr + 1; i <= nyr; i++){
	    M(ig)(i) = M(ig)(i - 1) * mfexp(log_m_devs(i));
	  }
	}

	// |---------------------------------------------------------------------------------|
	// | TOTAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	for(ig = 1; ig <= n_ags; ig++){
	  Z(ig) = M(ig) + F(ig);
	  S(ig) = mfexp(-Z(ig));
	}

	/*
	  Purpose:
	    This function initializes the numbers-at-age matrix in syr
	    based on log_rinit and log_init_rec_devs, the annual recruitment
	    based on log_rbar and log_rec_devs, and updates the number-at-age
	    over time based on the survival rate calculated in calcTotalMortality.

	  Notes:
	    - Aug 9, 2012.  Made a change here to initialize the numbers
	      at age in syr using the natural mortality rate at age in syr.
	      Prior to this the average (m_bar) rate was used, since this
	      has now changed with new projection control files.  Should only
	      affect models that were using time varying natural mortality.
	    - d_iscamCntrl(5) is a flag to start at unfished conditions, so set N(syr,sage) = ro
	*/
FUNCTION calcNumbersAtAge
	int ig, ih, kgear;
	N.initialize();
	bt.initialize();
	// Vulnerable biomass to all gears
	vbt.initialize();
	for(ig = 1; ig <= n_ags; ig++){
	  f = n_area(ig);
	  g = n_group(ig);
	  ih = pntr_ag(f, g);
	  dvar_vector lx(sage, nage);
	  dvar_vector tr(sage, nage);
	  lx(sage) = 1.0;
	  for(j = sage; j < nage; j++){
	    lx(j + 1) = lx(j) * exp(-M(ig)(syr)(j));
	  }
	  lx(nage) /= (1. - exp(-M(ig)(syr, nage)));
	  if(d_iscamCntrl(5)){
	    // initialize at unfished conditions.
	    tr = log(ro(g)) + log(lx);
	  }else{
	    tr(sage) = log_avgrec(ih) + log_rec_devs(ih)(syr);
	    tr(sage + 1, nage) = log_recinit(ih) + init_log_rec_devs(ih);
	    tr(sage + 1, nage) += log(lx(sage + 1, nage));
	  }
	  N(ig)(syr)(sage, nage) = 1.0 / nsex * mfexp(tr);
	  //LOG<<"tr\n"<<tr<<"\n";
	  //LOG<<"tr.indexmin() = "<<tr.indexmin()<<"\n";
	  //LOG<<"tr.indexmax() = "<<tr.indexmax()<<"\n";
	  log_rt(ih)(syr - nage + sage, syr) = tr.shift(syr - nage + sage);
	  //LOG<<"After tr\n"<<tr<<"\n";
	  //LOG<<"tr.indexmin() = "<<tr.indexmin()<<"\n";
	  //LOG<<"tr.indexmax() = "<<tr.indexmax()<<"\n";
	  //exit(1);
	  for(i = syr; i <= nyr; i++){
	    if(i > syr){
	      log_rt(ih)(i) = (log_avgrec(ih) + log_rec_devs(ih)(i));
 	      N(ig)(i, sage) = 1.0 / nsex * mfexp(log_rt(ih)(i));
	    }
	    N(ig)(i + 1)(sage + 1, nage) = ++elem_prod(N(ig)(i)(sage, nage - 1),
	                                               S(ig)(i)(sage, nage - 1));
	    N(ig)(i+1, nage) += N(ig)(i, nage) * S(ig)(i, nage);
	    // average biomass for group in year i
	    bt(g)(i) = sum(elem_prod(N(ig)(i), d3_wt_avg(ig)(i)));
	    // Vulnerable biomass to all gears
	    for(kgear = 1; kgear <= ngear; kgear++){
	      if(ig == 1){
	        vbt(g)(kgear)(i) = sum(elem_prod(elem_prod(N(ig)(i), d3_wt_avg(ig)(i)),
	                                     mfexp(log_sel(kgear)(ig)(i))));
	      }else{
	        vbt(g)(kgear)(i) += sum(elem_prod(elem_prod(N(ig)(i), d3_wt_avg(ig)(i)),
	                                     mfexp(log_sel(kgear)(ig)(i))));
	      }
	    }
	  }
	  N(ig)(nyr + 1, sage) = 1.0 / nsex * mfexp(log_avgrec(ih));
	  bt(g)(nyr + 1) = sum(elem_prod(N(ig)(nyr+1),
	                               d3_wt_avg(ig)(nyr+1)));
	  // Vulnerable biomass to all gears
	  for(kgear = 1; kgear <= ngear; kgear++){
	    if(ig == 1){
	      vbt(g)(kgear)(nyr + 1) = sum(elem_prod(elem_prod(N(ig)(nyr + 1),
	                                                       d3_wt_avg(ig)(nyr + 1)),
	                                             mfexp(log_sel(kgear)(ig)(nyr))));
	    }else{
	      vbt(g)(kgear)(nyr + 1) += sum(elem_prod(elem_prod(N(ig)(nyr + 1),
	                                                       d3_wt_avg(ig)(nyr + 1)),
	                                             mfexp(log_sel(kgear)(ig)(nyr))));
	    }
	  }
	}
	/*
	  Purpose:
	    This function calculates the predicted age-composition samples (A) for
	    both directed commercial fisheries and survey age-composition data. For
	    all years of data specified in the A matrix, calculated the predicted
	    proportions-at-age in the sampled catch-at-age.  If no catch-age data exist
	    for a particular year i, for gear k (i.e. no directed fishery or from a
	    survey sample process that does not have an appreciable F), the calculate
	    the predicted proportion based on log(N) + log_sel(group,gear,year)
	  Notes:
	    - Adapted from iSCAM 1.5.
	    - No longer using ragged arrays for gear, the ragged matrix is indexed by:
	      year gear area, group, sex | age columns
	    - For the residuals, note that each gear is weigthed by the conditional MLE
	      of the variance.
	*/
FUNCTION calcComposition
	int ii, ig, kk;
	ig = 0;
	dvar_vector va(sage, nage);
	dvar_vector fa(sage, nage);
	dvar_vector sa(sage, nage);
	dvar_vector za(sage, nage);
	dvar_vector ca(sage, nage);
	dvar_vector na(sage, nage);
	A_hat.initialize();
	for(kk = 1; kk <= nAgears; kk++){
	  for(ii = 1; ii <= n_A_nobs(kk); ii++){
	    i = d3_A(kk)(ii)(n_A_sage(kk) - 6);
	    k = d3_A(kk)(ii)(n_A_sage(kk) - 4);
	    f = d3_A(kk)(ii)(n_A_sage(kk) - 3);
	    g = d3_A(kk)(ii)(n_A_sage(kk) - 2);
	    h = d3_A(kk)(ii)(n_A_sage(kk) - 1);
	    // trap for retrospecitve analysis.
	    if(i < syr)
	      continue;
	    if(i > nyr)
	      continue;
	    if(h){
	      // age comps are sexed (h > 0)
	      ig = pntr_ags(f, g, h);
	      va = mfexp(log_sel(k)(ig)(i));
	      za = Z(ig)(i);
	      sa = S(ig)(i);
	      na = N(ig)(i);
	      if(ft(ig)(k)(i) == 0){
	        fa = va;
	      }else{
	        fa = ft(ig)(k)(i) * va;
	      }
	      ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);
	      //LOG<<"XXX\nZ("<<ig<<")("<<i<<") = "<<Z(ig)(i)<<"\n";
	      //LOG<<"S("<<ig<<")("<<i<<") = "<<S(ig)(i)<<"\n";
	      //LOG<<"N("<<ig<<")("<<i<<") = "<<N(ig)(i)<<"\n\n";
	      //  A_hat(kk)(ii) = ca(n_A_sage(kk),n_A_nage(kk));
	      //    +group if n_A_nage(kk) < nage
	      //
	      //  if(n_A_nage(kk) < nage){
	      //    A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
	      //  }
	    }else{
	      for(h = 1; h <= nsex; h++){
	        ig = pntr_ags(f,g,h);
	        va = mfexp(log_sel(k)(ig)(i));
	        za = Z(ig)(i);
	        sa = S(ig)(i);
	        na = N(ig)(i);
	        if(ft(ig)(k)(i) == 0){
	          fa = va;
	        }else{
	          fa = ft(ig)(k)(i) * va;
	        }
	        ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);
	        //  A_hat(kk)(ii) += ca(n_A_sage(kk),n_A_nage(kk));
	        //    +group if n_A_nage(kk) < nage
	        //  if(n_A_nage(kk) < nage){
	        //    A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
	        //  }
	      }
	    }
	    // This is the age-composition
	    if(n_ageFlag(kk)){
	      A_hat(kk)(ii) = ca(n_A_sage(kk), n_A_nage(kk));
	      if(n_A_nage(kk) < nage){
	        A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
	      }
	    }else{
	      /*
	        This the catch-at-length composition.
	        Pseudocode:
	        -make an ALK
	        -Ahat = ca * ALK
	      */
	      dvar_vector mu = d3_len_age(ig)(i);
	      dvar_vector sig = 0.1 * mu;
	      dvector x(n_A_sage(kk), n_A_nage(kk));
	      x.fill_seqadd(n_A_sage(kk), 1);
	      dvar_matrix alk = ALK(mu,sig,x);
	      A_hat(kk)(ii) = ca * alk;
	    }
	    A_hat(kk)(ii) /= sum( A_hat(kk)(ii));
	  }
	}

	/*
	  Purpose: This function calculates the total catch.
	  Dependencies: Must call calcCatchAtAge function first.
	  NOTES:
	  TODO list:
	    [ ] get rid of the obs_ct, ct, eta array structures, inefficient, better to use
	        a matrix, then cbind the predicted catch and residuals for report. (ie. an R
	        data.frame structure and use melt to ggplot for efficient plots.)
	*/
FUNCTION calcTotalCatch
	int ii, l, ig;
	double d_ct;
	ct.initialize();
	eta.initialize();
	dvar_vector fa(sage,nage);
	dvar_vector ca(sage,nage);
	dvar_vector sa(sage,nage);
	dvar_vector za(sage,nage);
	for(ii = 1; ii <= nCtNobs; ii++){
	  i = dCatchData(ii, 1); // year
	  k = dCatchData(ii, 2); // gear
	  f = dCatchData(ii, 3); // area
	  g = dCatchData(ii, 4); // group
	  h = dCatchData(ii, 5); // sex
	  l = dCatchData(ii, 6); // type
	  d_ct = dCatchData(ii, 7); // value
	  // | trap for retro year
	  if(i < syr || i > nyr)
	    continue;
	  switch(l){
	    case 1: // catch in weight
	      if(h){
	        ig = pntr_ags(f, g, h);
	        fa = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
	        za = Z(ig)(i);
	        sa = S(ig)(i);
	        ca = elem_prod(elem_prod(elem_div(fa, za), 1. - sa), N(ig)(i));
	        ct(ii) = ca * d3_wt_avg(ig)(i);
	      }else{
	        for(h = 1 ; h <= nsex; h++){
	          ig = pntr_ags(f, g, h);
	          fa = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
	          za = Z(ig)(i);
	          sa = S(ig)(i);
	          ca = elem_prod(elem_prod(elem_div(fa, za), 1. - sa), N(ig)(i));
	          ct(ii) += ca * d3_wt_avg(ig)(i);
	          //LOG<<"calcTotalCatch()\nft("<<ig<<")("<<k<<")("<<i<<") = "<<ft(ig)(k)(i)<<"\n";
	          //LOG<<"d3_wt_avg("<<ig<<")("<<i<<")\n"<<d3_wt_avg(ig)(i)<<"\n";
	          //LOG<<"ft("<<ig<<")("<<k<<")("<<i<<") = "<<ft(ig)(k)(i)<<"\n";
	          //LOG<<"fa\n"<<fa<<"\n";
	          //LOG<<"za\n"<<za<<"\n";
	          //LOG<<"sa\n"<<sa<<"\n";
	          //LOG<<"ca\n"<<ca<<"\n";
	          //LOG<<"ct("<<ii<<") = "<<ct(ii)<<"\n\n";
	        }
	      }
	      break;
	    case 2: // catch in numbers
	      if(h){
	        ig = pntr_ags(f, g, h);
	        fa = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
	        za = Z(ig)(i);
	        sa = S(ig)(i);
	        ca = elem_prod(elem_prod(elem_div(fa, za), 1. - sa), N(ig)(i));
	        ct(ii) = sum(ca);
	      }else{
	        for(h = 1; h <= nsex; h++){
	          ig = pntr_ags(f, g, h);
	          fa = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
	          za = Z(ig)(i);
	          sa = S(ig)(i);
	          ca = elem_prod(elem_prod(elem_div(fa, za), 1. - sa), N(ig)(i));
	          ct(ii) += sum(ca);
	        }
	      }
	      break;
	    case 3: // roe fisheries, special case
	      if(h){
	        ig = pntr_ags(f, g, h);
	        dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
	        ct(ii) = (1.-exp(-ft(ig)(k)(i))) * ssb;
	      }else{
	        for(h = 1; h <= nsex; h++){
	          ig = pntr_ags(f, g, h);
	          dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
	          ct(ii) += (1. - exp(-ft(ig)(k)(i))) * ssb;
	        }
	      }
	      break;
	  }
	  // | catch residual
	  eta(ii) = log(d_ct + TINY) - log(ct(ii) + TINY);
	}

	/*
	  Purpose: This function computes the mle for survey q, calculates the survey
	    residuals (epsilon).
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
	  TODO LIST:
	    [?] - add capability to accomodate priors for survey q\'s.
	    [ ] - verify q_prior=2 option for random walk in q.
	    [ ] - For sel_type==3, may need to reduce abundance by F on spawning biomass (herring)
	    [ ] - add capability to accommodate priors for survey catchability coefficients.
	*/
FUNCTION calcSurveyObservations
	int ii, kk, ig, nz;
	double di;
	dvariable ftmp;
	dvar_vector Na(sage,nage);
	dvar_vector va(sage,nage);
	dvar_vector sa(sage,nage);
	epsilon.initialize();
	std_epsilon.initialize();
	it_hat.initialize();
	for(kk = 1; kk <= nItNobs; kk++){
	  // Vulnerable number-at-age to survey.
	  dvar_matrix V(1,n_it_nobs(kk),sage,nage);
	  V.initialize();
	  nz = 0;
	  int iz = 1; // index for first year of data for prospective analysis.
	  for(ii = 1; ii <= n_it_nobs(kk); ii++){
	    i = d3_survey_data(kk)(ii)(1);
	    k = d3_survey_data(kk)(ii)(3);
	    f = d3_survey_data(kk)(ii)(4);
	    g = d3_survey_data(kk)(ii)(5);
	    h = d3_survey_data(kk)(ii)(6);
	    di = d3_survey_data(kk)(ii)(8);
	    // trap for retrospective nyr change
	    if(i < syr){
	      iz ++;
	      nz ++;
	      continue;
	    }
	    if(i > nyr)
	      continue;
	    nz ++; // counter for number of observations.
	    Na.initialize();
	    for(h = 1; h <= nsex; h++){
	      ig = pntr_ags(f,g,h);
	      va = mfexp(log_sel(k)(ig)(i));
	      sa = mfexp(-Z(ig)(i) * di);
	      Na = elem_prod(N(ig)(i),sa);
	      switch(n_survey_type(kk)){
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
	  }
	  dvector it = trans(d3_survey_data(kk))(2)(iz,nz);
	  dvector wt = log(trans(d3_survey_data(kk))(7)(iz,nz));
	  wt = wt / sum(wt);
	  dvar_vector t1 = rowsum(V);
	  dvar_vector zt = log(it) - log(t1(iz,nz));
	  dvariable zbar = sum(elem_prod(zt,wt));
	  q(kk) = mfexp(zbar);
	  // Survey residuals
	  epsilon(kk).sub(iz,nz) = zt - zbar;
	  std_epsilon(kk).sub(iz,nz) = (epsilon(kk).sub(iz,nz) - mean(epsilon(kk).sub(iz,nz))) /
	                                std_dev(epsilon(kk).sub(iz,nz));
	  it_hat(kk).sub(iz,nz) = q(kk) * t1(iz,nz);
	  // SPECIAL CASE: penalized random walk in q.
	  if(q_prior(kk) == 2){
	    epsilon(kk).initialize();
	    dvar_vector fd_zt = first_difference(zt);
	    dvariable zw_bar = sum(elem_prod(fd_zt,wt(iz,nz-1)));
	    epsilon(kk).sub(iz,nz-1) = fd_zt - zw_bar;
	    std_epsilon(kk).sub(iz,nz-1) = (epsilon(kk).sub(iz,nz-1) - mean(epsilon(kk).sub(iz,nz-1))) /
	                                    std_dev(epsilon(kk).sub(iz,nz-1));
	    qt(kk)(iz) = exp(zt(iz));
	    for(ii = iz + 1; ii <= nz; ii++){
	      qt(kk)(ii) = qt(kk)(ii-1) * exp(fd_zt(ii-1));
	    }
	    it_hat(kk).sub(iz,nz) = elem_prod(qt(kk)(iz,nz),t1(iz,nz));
	  }
	}

	/*
	  Purpose:
	    This function is used to derive the underlying stock-recruitment
	    relationship that is ultimately used in determining MSY-based reference
	    points. The objective of this function is to determine the appropriate
	    Ro, Bo and steepness values of either the Beverton-Holt or Ricker  Stock-
	  Recruitment Model:
	    Beverton-Holt Model
	    \f$ Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau) \f$
	    Ricker Model
	    \f$ Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau) \f$
	    The definition of a stock is based on group only. At this point, spawning biomass
	    from all areas for a given group is the assumed stock, and the resulting
	     recruitment is compared with the estimated recruits|group for all areas.
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
	int ig, ih;
	rt.initialize();
	sbt.initialize();
	delta.initialize();
	dvariable phib;
	dvector fa(sage,nage); //fecundity here
	dvar_vector stmp(sage,nage);
	dvar_vector ma(sage,nage);
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector lx(sage,nage);
	dvar_vector lw(sage,nage);
	for(g = 1; g <= ngroup; g++){
	  lx(sage) = 1.0;
	  lw(sage) = 1.0;
	  phib = 0;
	  for(f = 1; f <= narea; f++){
	    for(h = 1; h <= nsex; h++){
	      ig = pntr_ags(f,g,h);
	      // Step 1. average natural mortality rate at age.
	      // Step 2. calculate survivorship
	      for(j = sage; j <= nage; j++){
	        ma(j) = mean(trans(M(ig))(j));
	        fa(j) = mean(trans(d3_wt_mat(ig))(j));
	        if(j > sage){
	          lx(j) = lx(j - 1) * mfexp(-ma(j - 1));
	        }
	        lw(j) = lx(j) * mfexp(-ma(j) * d_iscamCntrl(13));
	      }
	      lx(nage) /= 1.0 - mfexp(-ma(nage));
	      lw(nage) /= 1.0 - mfexp(-ma(nage));
	      // Step 3. calculate average spawing biomass per recruit.
	      phib += (1.0 / nsex) * lw * fa;
	      LOG<<"FIT: sex = "<<h<<", phib:\n"<<phib<<"\n";

	      // Step 4. compute spawning biomass at time of spawning.
	      for(i = syr; i <= nyr; i++){
	        stmp = mfexp(-Z(ig)(i) * d_iscamCntrl(13));
	        sbt(g,i) += elem_prod(N(ig)(i), d3_wt_mat(ig)(i)) * stmp;
	      }
	      // Step 5. spawning biomass projection under natural mortality only.
	      stmp = mfexp(-M(ig)(nyr) * d_iscamCntrl(13));
	      sbt(g,nyr+1) += elem_prod(N(ig)(nyr+1),d3_wt_mat(ig)(i)) * stmp;
	    }
	    // | Estimated recruits
	    ih = pntr_ag(f,g);
	    rt(g) += mfexp(log_rt(ih)(syr+sage,nyr));
	  }
	  // Step 6. calculate stock recruitment parameters (so, beta, sbo);
	  so(g) = kappa(g) / phib;
	  sbo(g) = ro(g) * phib;
	  // Step 7. calculate predicted recruitment.
	  if(verbose){
	    //LOG<<"before g = "<<g<<", syr = "<<syr<<", sage = "<<sage<<", nyr = "<<nyr<<"\n";
	    //LOG<<"before sbt(g) = "<<sbt(g)<<"\n";
	  }
	  dvar_vector tmp_st = sbt(g)(syr,nyr-sage).shift(syr+sage);
	  if(verbose){
	    //LOG<<"after sbt(g) = "<<sbt(g)<<"\n\n";
	  }
	  switch(int(d_iscamCntrl(2))){
	    case 1: // Beverton Holt model
	      beta(g) = (kappa(g)-1.) / sbo(g);
	      tmp_rt = elem_div(so(g) * tmp_st, 1. + beta(g) * tmp_st);
	      break;
	    case 2: // Ricker model
	      beta(g) = log(kappa(g))/sbo(g);
	      tmp_rt = elem_prod(so(g) * tmp_st,exp(-beta(g) * tmp_st));
	      break;
	  }
	  // Step 8. // residuals in stock-recruitment curve with gamma_r = 0
	  delta(g) = log(rt(g)) - log(tmp_rt) + 0.5 * tau(g) * tau(g);
	  // Autocorrelation in recruitment residuals.
	  // if gamma_r > 0 then
	  if(active(gamma_r)){
	    int byr = syr + sage + 1;
	    delta(g)(byr,nyr) = log(rt(g)(byr,nyr))
	                          - (1.0-gamma_r) * log(tmp_rt(byr,nyr))
	                          - gamma_r * log(++rt(g)(byr-1,nyr-1))
	                          + 0.5 * tau(g) * tau(g);
	  }
	}

FUNCTION calcAnnualMeanWeight
	int ii, kk, ig, nz;
	double di;
	dvar_vector wNa(sage,nage);
	dvar_vector wva(sage,nage);
	dvar_vector wsa(sage,nage);
	wNa.initialize();
	wva.initialize();
	wsa.initialize();
	for(kk = 1; kk <= nMeanWt; kk++){ //loop through series with empirical annual mean weight data
	  dvar_matrix Vn(1,nMeanWtNobs(kk),sage,nage); // | Vulnerable number-at-age to gear
	  dvar_matrix Vb(1,nMeanWtNobs(kk),sage,nage); // | Vulnerable biomass-at-age to gear
	  Vn.initialize();
	  Vb.initialize();
	  nz = 0;
	  int iz = 1; // index for first year of data for prospective analysis.
	  for(ii = 1; ii <= nMeanWtNobs(kk); ii++){ //Loop through years
	    i = d3_mean_wt_data(kk)(ii)(1); //year
	    k = d3_mean_wt_data(kk)(ii)(3); //gear
	    f = d3_mean_wt_data(kk)(ii)(4); //area
	    g = d3_mean_wt_data(kk)(ii)(5); //group
	    h = d3_mean_wt_data(kk)(ii)(6); //sex
	    di = d3_mean_wt_data(kk)(ii)(7); //timing
	    // | trap for retrospective nyr change
	    if(i < syr){
	      iz++;
	      nz++;
	      continue;
	    }
	    if(i > nyr)
	      continue;
	    nz++; // counter for number of observations.
	    for(h = 1; h <= nsex; h++){
	      ig = pntr_ags(f,g,h);
	      wva = mfexp(log_sel(k)(ig)(i));
	      wsa = mfexp(-Z(ig)(i) * di); //accounts for survey timing
	      wNa = elem_prod(N(ig)(i),wsa);
	      Vn(ii) += elem_prod(wNa,wva); //adds sexes
	      Vb(ii) += elem_prod(elem_prod(wNa,wva),d3_wt_avg(ig)(i));
	    }
	    annual_mean_weight(kk)(ii) = sum(Vb(ii)) / sum(Vn(ii));
	    obs_annual_mean_weight(kk)(ii) = d3_mean_wt_data(kk)(ii)(2); //fill a matrix with observed annual mean weights - makes objective function calcs easier
	  }
	}

	/*
	  Purpose:  This function computes the objective function that ADMB will minimize.
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
	    [*] - Dec 20, 2010.  SJDM added prior to survey qs.
	    [ ] - q_prior is an ivector with current options of 0 & 1 & 2.
	      0 is a uniform density (ignored) and 1 is a normal
	      prior density applied to log(q), and 2 is a random walk in q.
	    [ ] - Allow for annual sig_c values in catch data likelihood.
	    [ ] - Increase dimensionality of sig and tau to ngroup.
	    [ ] - Correct likelihood for cases when rho > 0 (Schnute & Richards, 1995)
	    */
FUNCTION calcObjectiveFunction
	nlvec.initialize();
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR CATCH DATA
	// |---------------------------------------------------------------------------------|
	// | This likelihood changes between phases n-1 and n:
	// | Phase (n-1)  -> standard deviation in the catch based on user input d_iscamCntrl(3)
	// | Phase (n)    -> standard deviation in the catch based on user input d_iscamCntrl(4)
	// |
	double sig_c = d_iscamCntrl(3);
	if(last_phase()){
	  sig_c = d_iscamCntrl(4);
	}
	if(active(log_ft_pars)){
	  nlvec(1) = dnorm(eta, 0.0, sig_c);
	}

	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR RELATIVE ABUNDANCE INDICES
	// |---------------------------------------------------------------------------------|
	// | sig_it  -> vector of standard deviations based on relative wt for survey.
	for(k = 1; k <= nItNobs; k++){
	  ivector ig = it_grp(k);
	  dvar_vector sig_it(1,n_it_nobs(k));
	  for(i = 1; i <= n_it_nobs(k); i++){
	    sig_it(i) = sig(ig(i)) / it_wt(k,i);
	  }
	  nlvec(2,k) = dnorm(epsilon(k),sig_it);
	}

	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR AGE-COMPOSITION DATA
	// |---------------------------------------------------------------------------------|
	// | - Two options based on d_iscamCntrl(14):
	// | -  1 -> multivariate logistic using conditional MLE of the variance for weight.
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
	for(k = 1; k <= nAgears; k++){
	  if(n_A_nobs(k) > 0){
	    //int n_naa = 0; //retrospective counter
	    //int n_saa = 1; //prospective counter
	    int iyr;
	    dmatrix O(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
	    dvar_matrix P(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
	    dvar_matrix nu(n_saa(k),n_naa(k),n_A_sage(k),n_A_nage(k));
	    O.initialize();
	    P.initialize();
	    nu.initialize();
	    int ii = n_saa(k);
	    for(i = 1; i <= n_A_nobs(k); i++){
	      iyr = d3_A(k)(i)(n_A_sage(k) - 6); //index for year
	      if(iyr >= syr && iyr <= nyr){
	        O(ii) = d3_A_obs(k)(i).sub(n_A_sage(k),n_A_nage(k));
	        P(ii) = A_hat(k)(i).sub(n_A_sage(k),n_A_nage(k));
	        ii++;
	      }
	    }
	    // Choose form of the likelihood based on d_iscamCntrl(14) switch
	    // switch(int(d_iscamCntrl(14)))
	    logistic_normal cLN_Age(O,P,dMinP(k),dEps(k));
	    logistic_student_t cLST_Age(O,P,dMinP(k),dEps(k));
	    switch(int(nCompLikelihood(k))){
	      case 1:
	        nlvec(3,k) = dmvlogistic(O,P,nu,age_tau2(k),dMinP(k));
	        break;
	      case 2:
	        nlvec(3,k) = dmultinom(O,P,nu,age_tau2(k),dMinP(k));
	        break;
	      case 3:
	        if(!active(log_age_tau2(k))){
	          nlvec(3,k) = cLN_Age();
	        }else{
	          nlvec(3,k) = cLN_Age(exp(log_age_tau2(k)));
	        }
	        // Residual
	        if(last_phase()){
	          nu = cLN_Age.get_standardized_residuals();
	          age_tau2(k) = cLN_Age.get_sigma2();
	        }
	        break;
	      case 4:
	        //logistic_normal cLN_Age( O,P,dMinP(k),dEps(k) );
	        if(active(phi1(k)) && !active(phi2(k))){
	          //LOG<<"Running calcObjectiveFunction"<<'\n';
	          //LOG<<"        log_age_tau2: "<<log_age_tau2<<'\n';
	          //LOG<<"                   k: "<<k<<'\n';
	          //LOG<<"     log_age_tau2(k): "<<log_age_tau2(k)<<'\n';
	          //LOG<<"exp(log_age_tau2(k)): "<<exp(log_age_tau2(k))<<'\n';
	          //LOG<<"             phi2(k): "<<phi2(k)<<'\n';
	          //LOG<<" cLN_Age(expk,phi2k): "<<cLN_Age(exp(log_age_tau2(k)))<<'\n'<<'\n';
	          nlvec(3,k)   = cLN_Age(exp(log_age_tau2(k)),phi1(k));
	        }
	        if(active(phi1(k)) && active(phi2(k))){
	          nlvec(3,k) = cLN_Age(exp(log_age_tau2(k)),phi1(k),phi2(k));
	        }
	        // Residual
	        if(last_phase()){
	          nu = cLN_Age.get_standardized_residuals();
	          age_tau2(k) = cLN_Age.get_sigma2();
	        }
	        break;
	      case 5: // Logistic-normal with student-t
	        if(!active(log_degrees_of_freedom(k))){
	          nlvec(3,k) = cLST_Age();
	        }else{
	          nlvec(3,k) = cLST_Age(exp(log_degrees_of_freedom(k)));
	        }
	        // Residual
	        if(last_phase()){
	          nu = cLST_Age.get_standardized_residuals();
	          age_tau2(k) = cLST_Age.get_sigma2();
	        }
	        break;
	      case 6: // Multinomial with estimated effective sample size.
	        nlvec(3,k) = mult_likelihood(O,P,nu,log_degrees_of_freedom(k));
	        break;
	      case 7: // Multivariate-t
	        nlvec(3,k) = multivariate_t_likelihood(O,P,log_age_tau2(k), log_degrees_of_freedom(k), phi1(k),nu);
	        age_tau2(k) = exp(value(log_age_tau2(k)));
	        break;
	      case 8: // Dirichlet Multinomial
	        // Use Dirichlet Multinomial to estimate the predicted age matrices
	        //LOG<<"Gear "<<k<<"\nSample sizes\n"<<samp_sizes(k)<<"\n";exit(1);
	        for(int i = O.rowmin(); i <= O.rowmax(); i++){
	          temp_n = samp_sizes(k, i) * O(i);
	          nlvec(3,k) -= ddirmultinom(temp_n, P(i), log_phi(k));
	        }
	        if(last_phase()){
	          // Extract residuals.
	          int a = O.colmin();
	          int A = O.colmax();
	          int t = O.rowmin();
	          int T = O.rowmax();
	          dvariable Nsamp;
	          for(i = t; i<= T; i++){
	            int n = 0;
	            dvector oo = O(i) / sum(O(i));
	            dvar_vector pp = P(i) / sum(P(i));
	            // count number of observations greater than minp from control file (2% is a reasonable number)
	            for(int j = a; j <= A; j++)
	              if(oo(j) > dMinP(k))
	                n++;
	            ivector iiage(1,n);
	            dvector o1(1,n);
	            o1.initialize();
	            dvar_vector p1(1,n);
	            p1.initialize();
	            int kk = 1;
	            for(int j = a; j <= A; j++){
	              if(oo(j) <= dMinP(k)){
	                o1(kk) += oo(j);
	                p1(kk) += pp(j);
	              }else{
	                o1(kk) += oo(j);
	                p1(kk) += pp(j);
	                if(kk <= n)
	                  iiage(kk) = j; //ivector for the grouped residuals
	                if(kk < n)
	                  kk++;
	              }
	            }

	            // Variance for DM: Thorsen et. al 2017, Equation 8
	            // o1 is observed age comps data, p1 is estimated age comps
	            dvar_vector dm_variance = (elem_prod(o1, (1.0 - o1)) / samp_sizes(k,i)) *
	              ((samp_sizes(k,i) + exp(log_phi(k))) / (1.0 + exp(log_phi(k))));

	            dvar_vector pearson = elem_div(o1 - p1, sqrt(dm_variance));
	            for(int j = 1; j <= n; j++){
	              nu(i)(iiage(j)) = pearson(j);
	            }
	          }
	        }
	        break;
	    }
	    // Extract residuals.
	    for(i = n_saa(k); i <= n_naa(k); i++){
	      A_nu(k)(i)(n_A_sage(k),n_A_nage(k)) = nu(i);
	    }
	    ii = n_saa(k);
	    for(i = 1; i <= n_A_nobs(k); i++){
	      iyr = d3_A(k)(i)(n_A_sage(k) - 6); //index for year
	      if(iyr >= syr && iyr <= nyr){
	        A_nu(k)(i)(n_A_sage(k),n_A_nage(k)) = nu(ii++);
	        //LOG<<"iyr = "<<iyr<<", "<<d3_A(k)(i)<<"\n";
	        //LOG<<"iyr = "<<iyr<<", "<<A_nu(k)(i)(n_A_sage(k),n_A_nage(k))<<"\n";
	      }
	    }
	  }
	}

	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT LIKELIHOOD COMPONENT
	// |---------------------------------------------------------------------------------|
	// | - tau is the process error standard deviation.
	if(active(theta(1)) || active(theta(2))){
	  for(g = 1; g <= ngroup; g++){
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

	for(k = 1; k <= ngear; k++){
	  if(verbose){
	    //LOG<<"ACTIVE: k = "<<k<<", active(sel_par_f(k)) = "<<active(sel_par_f(k))<<", active(sel_par_m(k)) = "<<active(sel_par_m(k))<<"\n";
	  }
	  if(active(sel_par_f(k)) || active(sel_par_m(k))){
	    // If not using logistic selectivity then
	    if(isel_type(k) != 1 &&
	       isel_type(k) != 7 &&
	       isel_type(k) != 8 &&
	       isel_type(k) != 11){
	      for(ig = 1; ig <= n_ags; ig++){
	        for(i = syr; i <= nyr; i++){
	          //curvature in selectivity parameters
	          dvar_vector df2 = first_difference(first_difference(log_sel(k)(ig)(i)));
	          nlvec(5,k) += lambda_1(k) / (nage-sage+1) * df2 * df2;
	          //penalty for dome-shapeness
	          for(j = sage; j <= nage - 1; j++){
	            if(log_sel(k,ig,i,j) > log_sel(k,ig,i,j+1)){
	              nlvec(6,k) += lambda_2(k) * square(log_sel(k,ig,i,j) - log_sel(k,ig,i,j+1));
	            }
	          }
	        }
	      }
	      if(isel_type(k) == 4 || isel_type(k) == 5 || n_sel_blocks(k) > 1){
	        for(ig = 1; ig <= n_ags; ig++){
	          dvar_matrix trans_log_sel = trans(log_sel(k)(ig));
	          for(j = sage; j <= nage; j++){
	            dvar_vector df2 = first_difference(first_difference(trans_log_sel(j)));
	            nlvec(7,k) +=  lambda_3(k) / (nage - sage + 1) * norm2(df2);
	          }
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
	for(k = 1; k <= ngear; k++){
	  if(active(sel_par_f(k)) &&
	     isel_type(k) != 1 &&
	     isel_type(k) != 7 &&
	     isel_type(k) != 8 &&
	     isel_type(k) != 11){
	    dvariable s = 0;
	    if(isel_type(k) == 5){ //bicubic spline version ensure column mean = 0
	      dvar_matrix tmp = trans(sel_par_f(k));
	      for(j = 1; j <= tmp.rowmax(); j++){
	        s = mean(tmp(j));
	        lvec(1) += 10000.0 * s * s;
	      }
	    }
	    if(isel_type(k) == 2 ||
	       isel_type(k) == 3 ||
	       isel_type(k) == 4 ||
	       isel_type(k) == 12){
	      dvar_matrix tmp = sel_par_f(k);
	      for(j = 1; j <= tmp.rowmax(); j++){
	        s = mean(tmp(j));
	        lvec(1) += 10000.0 * s * s;
	      }
	    }
	  }
	  if(active(sel_par_m(k)) &&
	     isel_type(k) != 1 &&
	     isel_type(k) != 7 &&
	     isel_type(k) != 8 &&
	     isel_type(k) != 11){
	    dvariable s = 0;
	    if(isel_type(k) == 5){ //bicubic spline version ensure column mean = 0
	      dvar_matrix tmp = trans(sel_par_m(k));
	      for(j = 1; j <= tmp.rowmax(); j++){
	        s = mean(tmp(j));
	        lvec(1) += 10000.0 * s * s;
	      }
	    }
	    if(isel_type(k) == 2 ||
	       isel_type(k) == 3 ||
	       isel_type(k) == 4 ||
	       isel_type(k) == 12){
	      dvar_matrix tmp = sel_par_m(k);
	      for(j = 1; j <= tmp.rowmax(); j++){
	        s = mean(tmp(j));
	        lvec(1) += 10000.0 * s * s;
	      }
	    }
	  }
	}

	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR ANNUAL MEAN WEIGHT DATA
	// |---------------------------------------------------------------------------------|
	// | - sig_it     -> vector of standard deviations based on relative wt for survey.
	// |  init_3darray d3_mean_wt_data(1,nMeanWt,1,nMeanWtNobs,1,7)
	if(fitMeanWt){
	  for(k = 1; k <= nMeanWt; k++){
	    dvar_vector epsilon_wt = log(annual_mean_weight(k)) - log(obs_annual_mean_weight(k));
	    nlvec(8,k) = dnorm(epsilon_wt,weight_sig(k));
	  }
	}
	/* LOG<<"nlvec\n";
	for(i = 1; i <= nlvec.indexmax(); i++){
	  LOG<<"nlvec("<<i<<"):\n";
	  LOG<<nlvec(i)<<"\n";
	}
	LOG<<"End of nlvec "<<'\n';
	*/
	// |---------------------------------------------------------------------------------|
	// | PRIORS FOR LEADING PARAMETERS p(theta)
	// |---------------------------------------------------------------------------------|
	// | - theta_prior is a switch to determine which statistical distribution to use.
	// |
	dvariable ptmp;
	dvar_vector priors(1,npar);
	priors.initialize();
	//LOG<<"theta_control\n"<<theta_control<<"\n";
	//exit(1);
	for(i = 1; i <= npar; i++){
	  ptmp = 0;
	  for(j = 1; j <= ipar_vector(i); j++){
	    if(active(theta(i))){
	      switch(theta_prior(i)){
	        case 1: //normal
	          ptmp += dnorm(theta(i,j), theta_control(i,6), theta_control(i,7));
	          break;
	        case 2:
	          ptmp += dlnorm(theta(i,j), theta_control(i,6), theta_control(i,7));
	          break;
	        case 3: //beta distribution (0-1 scale)
	          double lb, ub;
	          lb = theta_lb(i);
	          ub = theta_ub(i);
	          ptmp += dbeta((theta(i,j)-lb) / (ub - lb), theta_control(i,6), theta_control(i,7));
	          break;
	        case 4: //gamma distribution
	          ptmp += dgamma(theta(i,j), theta_control(i,6), theta_control(i,7));
	          break;
	        default: //uniform density
	          ptmp += log(1. / (theta_control(i,3) - theta_control(i,2)));
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
	for(k = 1; k <= nits; k++){
	  if(q_prior(k) == 1){
	    qvec(k) = dnorm(log(q(k)), mu_log_q(k), sd_log_q(k));
	  }
	}

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
	if(last_phase()){
	  pvec(1) = dnorm(log_fbar, log(d_iscamCntrl(7)), d_iscamCntrl(9));
	  // | Penalty for log_rec_devs (large variance here)
	  for(g = 1; g <= n_ag; g++){
	    // http://api.admb-project.org/dnorm_8cpp_source.html
	    // Line 285 -> dvariable dnorm( const dvar_vector& residual, const double& std )
	    // LaTeX equation: n_{logrecdevs}(\frac{1}{2}log(2\pi) + log(\sigma)) +\frac{1}{2}\frac{\sum_{i=1}^{n_{logrecdevs}} logrecdevs[i]^2}{\sigma^2}
	    pvec(4) += dnorm(log_rec_devs(g), 2.0);
	    pvec(5) += dnorm(init_log_rec_devs(g), 2.0);
	    dvariable s = 0;
	    s = mean(log_rec_devs(g));
	    pvec(6) += 1.e5 * s * s;
	    s = mean(init_log_rec_devs(g));
	    pvec(7) += 1.e5 * s * s;
	  }
	}else{
	  pvec(1) = dnorm(log_fbar, log(d_iscamCntrl(7)), d_iscamCntrl(8));
	  //Penalty for log_rec_devs (CV ~ 0.0707) in early phases
	  for(g = 1; g <= n_ag; g++){
	    // http://api.admb-project.org/dvect3_8cpp_source.html
	    // Line 87 -> double norm2(const dvector& t1)
	    // Inside norm2() the multiplication operator '*' is overloaded for two dvectors:
	    // http://api.admb-project.org/dvector_8cpp_source.html
	    // Line 512 -> double operator*(const dvector& t1, const dvector& t2)
	    // LaTeX equation: 100\frac{\sum_{i=1}^{n_{logrecdevs}} logrecdevs[i]^2}
	    pvec(4) += 100.* norm2(log_rec_devs(g));
	    pvec(5) += 100.* norm2(init_log_rec_devs(g));
	    dvariable s = 0;
	    s = mean(log_rec_devs(g));
	    pvec(6) += 1.e5 * s * s;
	    s = mean(init_log_rec_devs(g));
	    pvec(7) += 1.e5 * s * s;
	  }
	}
	if(active(log_m_nodes)){
	  double std_mdev = d_iscamCntrl(11);
	  dvar_vector fd_mdevs = first_difference(log_m_devs);
	  pvec(2) = dnorm(fd_mdevs, std_mdev);
	  pvec(2) += 0.5 * norm2(log_m_nodes);
	}
	objfun  = sum(nlvec);
	objfun += sum(lvec);
	objfun += sum(priors);
	objfun += sum(pvec);
	objfun += sum(qvec);
	/*
	  LOG<<"priors\n"<<priors<<'\n';
	  LOG<<"pvec\n"<<pvec<<'\n';
	  LOG<<"qvec\n"<<qvec<<'\n';
	  LOG<<"objfun\n"<<objfun<<'\n';
	*/
	nf++;

	/*
	  Purpose: This function calculates the MSY-based reference points, and also loops
	    over values of F and F-multipliers and calculates the equilibrium yield
	    for each fishing gear.
	  NOTES:
	  PSEUDOCODE:
	    (1) : Construct array of selectivities (potentially sex based log_sel)
	    (2) : Construct arrays of d3_wt_avg and d3_wt_mat for reference years.
	    (3) : Come up with a reasonable guess for fmsy for each gear in nfleet.
	    (4) : Instantiate an Msy class object and get_fmsy.
	    (5) : Use Msy object to get reference points.
	  TODO list:
	    [ ] - allow user to specify which selectivity years are used in reference point
	          calculations. This should probably be done in the projection File Control.
	*/
FUNCTION void calcReferencePoints()
	if(d_iscamCntrl(13) || d_iscamCntrl(17)){
	  run_FRP(); // Code for testing reference point calcs
	}
	if(!d_iscamCntrl(13)){
	  int kk, ig;
	  // Matrix of selectivities for directed fisheries.
	  // log_sel(gear)(n_ags)(year)(age)
	  // ensure dAllocation sums to 1.
	  dvector d_ak(1,nfleet);
	  d3_array d_V(1,n_ags,1,nfleet,sage,nage);
	  dvar3_array dvar_V(1,n_ags,1,nfleet,sage,nage);
	  for(k = 1;k <= nfleet;k++){
	    kk = nFleetIndex(k);
	    d_ak(k) = dAllocation(kk);
	    if(verbose){
	      //LOG<<"kk = nFleetIndex(k) = "<<nFleetIndex(k)<<"\n";
	      //LOG<<"dAllocation(kk) = "<<dAllocation(kk)<<"\n\n";
	    }
	    for(ig = 1;ig <= n_ags;ig++){
	      d_V(ig)(k) = value(exp(log_sel(kk)(ig)(nyr)));
	      dvar_V(ig)(k) = exp(log_sel(kk)(ig)(nyr));
	    }
	  }
	  if(verbose){
	    //LOG<<"d_ak\n"<<d_ak<<"\n";
	  }
	  d_ak /= sum(d_ak);
	  // Average weight and mature spawning biomass for reference years
	  // dWt_bar(1,n_ags,sage,nage)
	  dmatrix fa_bar(1,n_ags,sage,nage);
	  dmatrix M_bar(1,n_ags,sage,nage);
	  for(ig = 1; ig <= n_ags; ig++){
	    fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
	    M_bar(ig)  = colsum(value(M(ig).sub(pf_cntrl(1),pf_cntrl(2))));
	    M_bar(ig) /= pf_cntrl(2) - pf_cntrl(1) + 1;
	  }
	  // Initial guess for fmsy for each fleet
	  // set fmsy = 2/3 of M divided by the number of fleets
	  fmsy.initialize();
	  fall.initialize();
	  msy.initialize();
	  bmsy.initialize();
	  dvar_vector dftry(1,nfleet);
	  dftry = 0.6/nfleet * mean(M_bar);
	  // Instantiate msy class for each stock
	  if(verbose & last_phase()){
	    for(g = 1; g <= ngroup; g++){
	      // Check that average weights are the same in calcReferencePoints as in the slow_MSY code
	      //double d_rho = d_iscamCntrl(13);
	      dvector d_mbar = M_bar(g);
	      dvector d_wa = dWt_bar(g);
	      dvector d_fa = fa_bar(g);
	      //LOG<<"g = "<<g<<"\n";
	      //LOG<<"Weight-at-age - dWt_bar(g)\n"<<dWt_bar(g)<<"\n";
	      //LOG<<"Fecundity-at-age - fa_bar(g)\n"<<fa_bar(g)<<"\n";
	      //LOG<<"One sex or Female Mean Natural mortality - M_bar(1)\n"<<M_bar(1)<<"\n";
	      if(nsex == 2){
	        //LOG<<"Male Natural mortality - M(2)\n"<<M(2)<<"\n";
	      }
	      //LOG<<"nfleet = "<<nfleet<<"\n";
	      //LOG<<"mean(M_bar) = "<<mean(M_bar)<<"\n";
	      //LOG<<"dftry = 0.6/nfleet * mean(M_bar) = "<<dftry<<"\n";
	      //LOG<<"Selectivity - dvar_V:\n"<<dvar_V<<"\n";
	      //LOG<<"ro(g) = "<<ro(g)<<"\n";
	      //LOG<<"steepness(g) = "<<steepness(g)<<"\n";
	      //LOG<<"d_rho = "<<d_rho<<"\n\n";
	    }
	  }
	  for(ig = 1; ig <= n_ags; ig++){
	    fa_bar(ig) = elem_prod(dWt_bar(ig), ma(ig));
	    M_bar(ig) = colsum(value(M(ig).sub(pf_cntrl(1), pf_cntrl(2))));
	    M_bar(ig) /= pf_cntrl(2) - pf_cntrl(1) + 1;
	  }
	  for(g = 1; g <= ngroup; g++){
	    double d_ro = value(ro(g));
	    double d_h = value(steepness(g));
	    double d_rho = d_iscamCntrl(13);
	    Msy c_msy(d_ro, d_h, M_bar, d_rho, dWt_bar, fa_bar, &d_V);
	    fmsy(g) = 0.1;
	    c_msy.get_fmsy(fmsy(g), d_ak);
	    bmsy(g) = c_msy.getBmsy();
	    msy(g) = c_msy.getMsy();
	    bo = c_msy.getBo();
	  }
	}

	/*
	  - This is a simple test routine for comparing the MSY class output to the
	  - MSF.xlsx spreadsheet that was used to develop the multiple fleet msy
	  - method. Its a permanent feature of the iscam code for testing.
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
	//LOG<<"Initial Fe "<<dftry<<'\n';
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSY(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dfmsy = c_MSY.getFmsy(dftry);
	//LOG<<"Fmsy = "<<dfmsy<<'\n';
	dvector ak(1,2);
	ak = 0.3;
	ak(2) = 1-ak(1);
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSYk(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dkmsy = c_MSYk.getFmsy(dftry,ak);
	//LOG<<"Fmsy_k ="<<dkmsy<<'\n';
	c_MSYk.print();
	dvector akmsy = c_MSYk.getFmsy(dftry);
	c_MSYk.print();
	exit(1);

	/*
	  Purpose: This function returns a modified selectivity vector (va) based on
	    the assumption that age-based selectivity will operate on the principle
	    of ideal free distribution.
	  Arguments:
	    \param  va -> age-specific vulnerability
	    \param  ba -> age-specific biomass (relative abundance is fine)
	  NOTES:
	  TODO list:
	    [ ]
	*/
FUNCTION dvector ifdSelex(const dvector& va, const dvector& ba, const double& mpow)
	dvector pa(sage,nage);
	pa = (elem_prod(va,pow(ba,mpow)));
	pa = pa/sum(pa);
	pa = exp( log(pa) - log(mean(pa)) );
	return (pa);

REPORT_SECTION
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
	report<<"init_log_rec_devs\n"<<init_log_rec_devs<<"\n";
	report<<"log_rec_devs\n"<<log_rec_devs<<"\n";;
	report<<"like_catch\n"<<nlvec(1)<<"\n";
	report<<"like_survey_index\n"<<nlvec(2)<<"\n";
	report<<"like_age_comps\n"<<nlvec(3)<<"\n";
	report<<"like_stock_recruit\n"<<nlvec(4)<<"\n";
	report<<"like_fishery_sel_curvature\n"<<nlvec(5)<<"\n";
	report<<"like_fishery_sel_dome_shapedness\n"<<nlvec(6)<<"\n";
	report<<"like_fishery_sel_first_differences\n"<<nlvec(7)<<"\n";
	report<<"like_annual_mean_weights\n"<<nlvec(8)<<"\n";
	REPORT(ro);
	dvector rbar = value(exp(log_avgrec));
	REPORT(rbar);
	dvector rinit = value(exp(log_recinit));
	REPORT(rinit);
	REPORT(sbo);
	REPORT(kappa);
	dvector steepness = value(theta(2));
	REPORT(steepness);
	REPORT(m);
	// double tau = value(sqrt(1.-rho)*varphi);
	// double sig = value(sqrt(rho)*varphi);
	report<<"rho\n"<<theta(6 + nsex - 1)<<'\n';
	report<<"vartheta\n"<<theta(7 + nsex - 1)<<'\n';
	REPORT(varphi);
	report<<"log_phi\n";
	for(k = 1; k <= nAgears; k++){
	  report<<log_phi(k)<<"\n";
	}
	report<<"dm_sample_sizes\n";
	for(k = 1; k <= nAgears; k++){
	  report<<samp_sizes(k)<<"\n";
	}
	report<<"Index SDs (sig_it)\n";
	for(k = 1; k <= nItNobs; k++){
	  ivector ig = it_grp(k);
	  //dvar_vector sig_it(1,n_it_nobs(k));
	  report<<"Gear "<<k<<":\n";
	  for(i = 1; i <= n_it_nobs(k); i++){
	    report<<sig(ig(i)) / it_wt(k,i)<<"\n";
	  }
	}

	REPORT(tau);
	REPORT(sig);
	REPORT(age_tau2);
	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS & AGE-SCHEDULE INFORMATION ON GROWTH AND MATURITY
	// |---------------------------------------------------------------------------------|
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
	REPORT(age);
	REPORT(la);
	REPORT(wa);
	REPORT(ma);
	// |---------------------------------------------------------------------------------|
	// | OBSERVED AND PREDICTED DATA AND RESIDUALS
	// |---------------------------------------------------------------------------------|
	// | Catch data
	// | Survey data
	// | Age composition data
	// | Empirical weight-at-age data
	REPORT(dCatchData);
	report<<"ct"<<'\n';
	int iter = 1;
        for(int fleet = 1; fleet <= nfleet; fleet++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(fleet != 1 && yr == syr){
	      report<<"\n";
	    }else if(iter > 1){
	      report<<" ";
	    }
	    report<<ct(iter);
	    iter++;
	  }
	}
		report<<"\n";

	REPORT(eta);
	REPORT(q);
	REPORT(qt);
	REPORT(d3_survey_data);
	REPORT(it_hat);
	REPORT(it_wt);
	REPORT(epsilon);
	if(n_A_nobs(1) > 0){
	  REPORT(n_A_sage);
	  REPORT(n_A_nage);
	  for(k = 1; k <= nAgears; k++){
	    report<<"d3_A_gear_"<<column(d3_A(k), -3)(1)<<"\n";
	    report<<d3_A(k)<<"\n";
	  }
	  // A_hat and A_nu have the same structure as the observed age comps, d3_A
	  // so prepend the headers from d3_A to each row of output for them.
	  int year, gear, area, group, sex;
	  for(k = 1; k <= nAgears; k++){
	    report<<"A_hat_gear_"<<column(d3_A(k), -3)(1)<<"\n";
	    for(int row = A_hat(k).rowmin(); row <= A_hat(k).rowmax(); row++){
	      year = d3_A(k, row)(-5);
	      gear = d3_A(k, row)(-3);
	      area = d3_A(k, row)(-2);
	      group = d3_A(k, row)(-1);
	      sex = d3_A(k, row)(0);
	      report<<" "<<year<<" "<<gear<<" "<<area<<" "<<group<<" "<<sex<<A_hat(k, row)<<"\n";
	    }
	  }
	  for(k = 1; k <= nAgears; k++){
	    report<<"A_nu_gear_"<<column(d3_A(k), -3)(1)<<"\n";
	    for(int row = A_nu(k).rowmin(); row <= A_nu(k).rowmax(); row++){
	      year = d3_A(k, row)(-5);
	      gear = d3_A(k, row)(-3);
	      area = d3_A(k, row)(-2);
	      group = d3_A(k, row)(-1);
	      sex = d3_A(k, row)(0);
	      report<<" "<<year<<" "<<gear<<" "<<area<<" "<<group<<" "<<sex<<A_nu(k, row)<<"\n";
	    }
	  }
	  report<<"Neff_multinomial\n";
	  dvector nscaler(1,nAgears);
	  nscaler.initialize();
	  int naa;
	  int iyr;
	  for(k = 1; k <= nAgears; k++){
	    if(int(nCompLikelihood(k))){
	      naa = 0;
	      //retrospective counter
	      for(i = 1; i <= n_A_nobs(k); i++){
	        iyr = d3_A(k)(i)(n_A_sage(k) - 6); //index for year
	        if(iyr<=nyr)
	          naa++;
	        else
	          continue;
	      }
	      dmatrix O = trans(trans(d3_A_obs(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
	      dvar_matrix P = trans(trans(A_hat(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
	      for(j = 1; j<= naa; j++){
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
	REPORT(dWt_bar);
	REPORT(d3_wt_mat);
	REPORT(d3_wt_dev);
	report<<"d3_wt_avg"<<'\n';
	for(int ig = 1; ig <= n_ags; ig++){
	  f = n_area(ig);
	  g = n_group(ig);
	  h = n_sex(ig);
	  for(i = syr; i <= nyr; i++){
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
	report<<"sel_par_f"<<'\n';
	for(k = 1; k <= ngear; k++){
	  for(j = 1; j <= jsel_npar(k); j++){
	    report<<k<<"\t"<<j<<"\t"<<exp(sel_par_f(k)(j))<<'\n';
	  }
	}
	report<<"sel_par_m"<<'\n';
	for(k = 1; k <= ngear; k++){
	  for(j = 1; j <= jsel_npar(k); j++){
	    report<<k<<"\t"<<j<<"\t"<<exp(sel_par_m(k)(j))<<'\n';
	  }
	}
	report<<"log_sel"<<'\n';
	for(k = 1; k <= ngear; k++){
	  for(int ig = 1; ig <= n_ags; ig++){
	    for(i = syr; i <= nyr; i++){
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
	for(int ig = 1; ig <= n_ags; ig++ ){
	  report<<ft(ig)<<'\n';
	}
	report<<"ut"<<'\n';
	for(int ig = 1; ig <= n_ags; ig++ ){
	  report<<1.0-exp(-ft(ig))<<'\n';
	}
	REPORT(M);
	REPORT(F);
	REPORT(Z);

	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT
	// |---------------------------------------------------------------------------------|
	int rectype = int(d_iscamCntrl(2));
	REPORT(rectype);
	REPORT(so);
	REPORT(beta);
	REPORT(sbt);
	REPORT(bt);
	REPORT(rt);
	REPORT(delta);
	report<<"vbt"<<'\n';
	for(k = 1; k <= ngear; k++){
	  for(int ig = 1; ig <= ngroup; ig++){
	    for(i = syr; i <= nyr + 1; i++){
	      report<<k<<"\t"<<ig<<"\t"<<i<<"\t"<<vbt(ig)(k)(i)<<'\n';
	    }
	  }
	}
	dmatrix rep_rt = value(exp(trans(trans(log_rt).sub(syr,nyr))));
	for(int ig = 1; ig <= n_ag; ig++){
	  rep_rt(ig)(syr) = value(exp( log_rt(ig)(syr-nage+sage)));
	}
	REPORT(rep_rt);
	// |---------------------------------------------------------------------------------|
	// | ABUNDANCE IN NUMBERS
	// |---------------------------------------------------------------------------------|
	REPORT(N);
	// |---------------------------------------------------------------------------------|
	// | ANNUAL MEAN WEIGHT DATA
	// |---------------------------------------------------------------------------------|
	REPORT(obs_annual_mean_weight);
	REPORT(annual_mean_weight);
	// |---------------------------------------------------------------------------------|
	// | MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	if(last_phase()){
	  if(verbose){
	    //LOG<<"\n\nCalculating MSY-based reference points\n";
	  }
	  calcReferencePoints();
	  if(verbose){
	    //LOG<<"Finished calcReferencePoints\n\n";
	  }
	  REPORT(bo);
	  REPORT(fmsy);
	  REPORT(msy);
	  REPORT(bmsy);
	  // REPORT(Umsy);
	  //LOG<<"Running Projections\n";
	  int ii;
	  for(ii = 1; ii <= n_tac; ii++){
	    //LOG<<"TAC "<<ii<<" = "<<tac(ii)<<'\n';
	    projection_model(tac(ii));
	  }
	  LOG<<" ______________________________ \n";
	  LOG<<"|    END OF REPORT SECTION     |\n";
	  LOG<<"|______________________________|\n";
	}
	// |---------------------------------------------------------------------------------|
	// | OUTPUT FOR OPERATING MODEL
	// |---------------------------------------------------------------------------------|
	// | Move to final section?
	if(last_phase()){
	  ofstream ofs("iSCAM.res");
	  ofs<<"# Bo\n"<<bo<<'\n';
	  ofs<<"# Fmsy\n"<<fmsy<<'\n';
	  ofs<<"# MSY\n"<<msy<<'\n';
	  ofs<<"# Bmsy\n"<<bmsy<<'\n';
	  ofs<<"# Sbt\n";
	  for(g = 1; g <= ngroup; g++){
	    ofs<<sbt(g)(nyr+1)<<"\t";
	  }
	  ofs<<'\n';
	  // projected biomass
	  // The total biomass for each stock
	  ofs<<"# Total biomass\n";
	  for(g = 1; g <= ngroup; g++){
	    ofs<<bt(g)(nyr+1)<<"\t";
	  }
	  ofs<<'\n';
	  ofs<<"# Numbers-at-age\n";
	  for(int ig = 1; ig <= n_ags; ig++){
	    ofs<<N(ig)(nyr+1)<<'\n';
	  }
	  ofs<<"# Weight-at-age\n";
	  for(int ig = 1; ig <= n_ags; ig++){
	    ofs<<d3_wt_avg(ig)(nyr+1)<<'\n';
	  }
	  ofs<<"# Natural mortality-at-age\n";
	  for(int ig = 1; ig <= n_ags; ig++ ){
	    ofs<<M(ig)(nyr)<<'\n';
	  }
	  // 4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
	  ofs<<"# log_selectivity\n";
	  for(int k = 1; k <= ngear; k++){
	    for(int ig = 1; ig <= n_ags; ig++){
	      ofs<<log_sel(k)(ig)(nyr)<<'\n';
	    }
	  }
	}
	/*
	  Stock status info
	  Bstatus = sbt / bmsy;
	  Fstatus = ft / fmsy;
	  if(bmsy > 0){
	    dvector Bstatus = value(sbt / bmsy);
	    REPORT(Bstatus);
	  }
	  dmatrix Fstatus(1,ngear,syr,nyr);
	  Fstatus.initialize();
	  for(k = 1; k <= nfleet; k++){
	    if(fmsy(k) >0){
	      j = nFleetIndex(k);
	      Fstatus(j) = value(ft(j) / fmsy(k));
	    }
	  }
	  REPORT(Fstatus);
	*/

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
FUNCTION decision_table
	int i;
	for(i = 1; i <= n_tac; i++){
	  LOG<<i<<" "<<tac<<'\n';
	  projection_model(tac(i));
	}
	LOG<<"Ok to here"<<'\n';

FUNCTION mcmc_output
	int iter;
	static int post_num = 0;
	post_num++;

	if(nf == 1){
	  LOG<<"In mcmc_output()\n";

	  // Open the files and write the headers
	  ofstream ofs("iscam_mcmc.csv");
	  // The structure for these objects can be found at roughly lines 924 and 1409.
	  // they are set up as vector_vectors to increase dimensionality
	  // Assume only one group and area for simplicity

	  ofs<<"ro";
	  ofs<<","<<"h";
	  for(int sex = 1; sex <= n_gs; sex++){
	    ofs<<","<<"m_sex"<<sex;
	  }
	  ofs<<","<<"rbar";
	  ofs<<","<<"rinit";
	  ofs<<","<<"rho";
	  ofs<<","<<"vartheta";
	  ofs<<","<<"bo";
	  ofs<<","<<"sbo";
	  if(!d_iscamCntrl(13) || (d_iscamCntrl(13) && !d_iscamCntrl(20))){
	    ofs<<","<<"bmsy";
	    for(int fleet = 1; fleet <= nfleet; fleet++){
	      ofs<<","<<"msy_fleet"<<fleet;
	      ofs<<","<<"fmsy_fleet"<<fleet;
	      ofs<<","<<"umsy_fleet"<<fleet;
	    }
	    for(int gear = 1; gear <= nItNobs; gear++){
	      ofs<<","<<"q_gear"<<gear;
	    }
	    ofs<<","<<"SSB";
	    for(k = 1; k <= ngear; k++){
	      for(int sel_blk = 1; sel_blk <= n_sel_blocks(k); sel_blk++){
	        ofs<<","<<"sel_age50_female_gear"<<k<<"_block"<<sel_blk;
	        ofs<<","<<"sel_sd50_female_gear"<<k<<"_block"<<sel_blk;
	        ofs<<","<<"sel_age50_male_gear"<<k<<"_block"<<sel_blk;
	        ofs<<","<<"sel_sd50_male_gear"<<k<<"_block"<<sel_blk;
	      }
	    }
	    ofs<<","<<"f";
	    ofs<<'\n';
	    ofstream of1("iscam_sbt_mcmc.csv");
	    for(int yr = syr; yr <= nyr + 1; yr++){
	      if(yr == syr){
	        of1<<"sbt"<<"_"<<yr;
	      }else{
	        of1<<",sbt"<<"_"<<yr;
	      }
	    }
	    of1<<'\n';
	    LOG<<""<<"\n";
	    ofstream of_ct("iscam_ct_mcmc.csv");
	    iter = 1;
	    for(int fleet = 1; fleet <= nfleet; fleet++){
	      for(int yr = syr; yr <= nyr; yr++){
	        if(iter > 1){
	          of_ct<<",";
	        }
	        of_ct<<"ct_fleet"<<fleet<<"_"<<yr;
		iter++;
	      }
	    }
	    of_ct<<"\n";
	    ofstream of2("iscam_rt_mcmc.csv");
	    for(int yr = syr + sage; yr <= nyr; yr++){
	      if(yr == syr + sage){
	        of2<<"rt"<<"_"<<yr;
	      }else{
	        of2<<",rt"<<"_"<<yr;
	      }
	    }
	    of2<<'\n';
	    ofstream of3("iscam_ft_mcmc.csv");
	    iter = 1;
	    for(int fleet = 1; fleet <= nfleet; fleet++){
	      for(int yr = syr; yr <= nyr; yr++){
	        if(iter == 1){
	          of3<<"ft_fleet"<<fleet<<"_"<<yr;
	        }else{
	          of3<<",ft_fleet"<<fleet<<"_"<<yr;
	        }
	        iter++;
	      }
	    }
	    of3<<'\n';
	    ofstream of4("iscam_rdev_mcmc.csv");
	    iter = 1;
	    for(int yr = syr + sage; yr <= nyr; yr++){
	      if(iter == 1){
	        of4<<"rdev_"<<yr;
	      }else{
	        of4<<",rdev_"<<yr;
	      }
	      iter++;
	    }
	    of4<<'\n';
	    ofstream of5("iscam_vbt_mcmc.csv");
	    iter = 1;
	    for(int fleet = 1; fleet <= nfleet; fleet++){
	      for(int yr = syr; yr <= nyr; yr++){
	        if(iter == 1){
	          of5<<"vbt_fleet"<<fleet<<"_"<<yr;
	        }else{
	          of5<<",vbt_fleet"<<fleet<<"_"<<yr;
	        }
	        iter++;
	      }
	    }
	    of5<<'\n';
	    ofstream of6("iscam_ut_mcmc.csv");
	    iter = 1;
	    iter = 1;
	    for(int fleet = 1; fleet <= nfleet; fleet++){
	      for(int yr = syr; yr <= nyr + 1; yr++){
	        if(iter == 1){
	          of6<<"ut_fleet"<<fleet<<"_"<<yr;
	        }else{
	          of6<<",ut_fleet"<<fleet<<"_"<<yr;
	        }
	        iter++;
	      }
	    }
	    of6<<'\n';
	    ofstream of7("iscam_m_mcmc.csv");
	    iter = 1;
	    for(int sex = 1; sex <= n_gs; sex++){
	      for(int yr = syr ; yr <= nyr; yr++){
	        if(iter == 1){
	          of7<<"m_sex"<<sex<<"_"<<yr;
	        }else{
	          of7<<",m_sex"<<sex<<"_"<<yr;
	        }
	        iter++;
	      }
	    }
	    of7<<'\n';
	    ofstream of8("iscam_index_fits_mcmc.csv");
	    iter = 1;
	    for(int kk = 1; kk <= nItNobs; kk++){
	      for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	        if(iter == 1){
	          of8<<"gear"<<kk<<"_yr"<<yr;
	        }else{
	          of8<<",gear"<<kk<<"_yr"<<yr;
	        }
	        iter++;
	      }
	    }
	    of8<<'\n';
	    ofstream of9("iscam_index_residuals_mcmc.csv");
	    iter = 1;
	    for(int kk = 1; kk <= nItNobs; kk++){
	      for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	        if(iter == 1){
	          of9<<"gear"<<kk<<"_yr"<<yr;
	        }else{
	          of9<<",gear"<<kk<<"_yr"<<yr;
	        }
	        iter++;
	      }
	    }
	    of9<<'\n';
	    ofstream of99("iscam_index_standardized_residuals_mcmc.csv");
	    iter = 1;
	    for(int kk = 1; kk <= nItNobs; kk++){
	      for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	        if(iter == 1){
	          of99<<"gear"<<kk<<"_yr"<<yr;
	        }else{
	          of99<<",gear"<<kk<<"_yr"<<yr;
	        }
	        iter++;
	      }
	    }
	    of99<<'\n';
	    ofstream of10("iscam_age_fits_mcmc.csv");
	    of10<<"gear,posterior,year,sex,";
	    for(int a = sage; a < nage; a++){
	      of10<<a<<",";
	    }
	    of10<<nage<<"\n";
	    ofstream of11("iscam_age_residuals_mcmc.csv");
	    of11<<"gear,posterior,year,sex,";
	    for(int a = sage; a < nage; a++){
	      of11<<a<<",";
	    }
	    of11<<nage<<"\n";
	    ofstream of12("iscam_selectivity_mcmc.csv");
	    of12<<"gear,posterior,block,start_year,end_year,sex,a_hat,g_hat\n";
	  }
	}
	// ---------------------------------------------------------------------
	// Headers done, output values
	// ---------------------------------------------------------------------
	// Leading parameters & reference points
	if(d_iscamCntrl(13)){
	  if(d_iscamCntrl(20)){
	    // Only report B0 - calc_bo() is in slowmsy.cpp
	    double B0;
	    calc_bo(B0, sage, nage, M, dWt_bar, ma, ro, d_iscamCntrl, pf_cntrl);
	    bo = B0;
	  }else{
	    calcReferencePoints();
	  }
	}else{
	  calcReferencePoints();
	}

	// Append the values to the files
	//these should be named parameters
	ofstream ofs("iscam_mcmc.csv",ios::app);
	ofs<<exp(theta(1)(1));
	ofs<<","<<theta(2)(1);
	ofs<<","<<exp(theta(3)(1));
	if(nsex == 2){
	  ofs<<","<<exp(theta(4)(1));
	}
	ofs<<","<<exp(log_avgrec(1));
	ofs<<","<<exp(theta(5 + nsex - 1)(1));
	ofs<<","<<theta(6 + nsex - 1)(1);
	ofs<<","<<theta(7 + nsex - 1)(1);
	ofs<<","<<bo;
	ofs<<","<<sbo;
	if(!d_iscamCntrl(13) || (d_iscamCntrl(13) && !d_iscamCntrl(20))){
	  ofs<<","<<bmsy;
	  for(int fleet = 1; fleet <= nfleet; fleet++){
	    ofs<<","<<msy(1,fleet);
	    ofs<<","<<fmsy(1,fleet);
	    ofs<<","<<1.0 - exp(-fmsy(1,fleet));
	  }
	}
	for(int gear = 1; gear <= nItNobs; gear++){
	  ofs<<","<<q(gear);
	}
	ofs<<","<<sbt(1)(nyr);
	for(k = 1; k <= ngear; k++){
	  for(int sel_blk = 1; sel_blk <= n_sel_blocks(k); sel_blk++){
	    ofs<<","<<exp(sel_par_f(k)(sel_blk)(1));
	    ofs<<","<<exp(sel_par_f(k)(sel_blk)(2));
	    ofs<<","<<exp(sel_par_m(k)(sel_blk)(1));
	    ofs<<","<<exp(sel_par_m(k)(sel_blk)(2));
	  }
	}
	ofs<<","<<objfun;
	ofs<<'\n';

	// output spawning stock biomass
	ofstream of1("iscam_sbt_mcmc.csv", ios::app);
	for(int yr = syr; yr <= nyr + 1; yr++){
	  if(yr == syr){
	    of1<<sbt(1)(yr);
	  }else{
	    of1<<","<<sbt(1)(yr);
	  }
	}
	of1<<'\n';

	// output spawning stock biomass
	ofstream of_ct("iscam_ct_mcmc.csv", ios::app);
	iter = 1;
        for(int fleet = 1; fleet <= nfleet; fleet++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(iter > 1){
	      of_ct<<",";
	    }
	    of_ct<<ct(iter);
	    iter++;
	  }
	}
	of_ct<<"\n";

	// output age-1 recruits
	ofstream of2("iscam_rt_mcmc.csv", ios::app);
	for(int yr = syr + sage; yr <= nyr; yr++){
	  if(yr == syr + sage){
	    of2<<rt(1)(yr);
	  }else{
	    of2<<","<<rt(1)(yr);
	  }
	}
	of2<<'\n';
	// output fishing mortality
	ofstream of3("iscam_ft_mcmc.csv", ios::app);
	iter = 1;
	for(int fleet = 1; fleet <= nfleet; fleet++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(iter == 1){
	      of3<<ft(1)(fleet)(yr);
	    }else{
	      of3<<","<<ft(1)(fleet)(yr);
	    }
	    iter++;
	  }
	}
	of3<<'\n';
	// output recruitment deviations
	// This is what the declaration of log_dev_recs looks like:
	// init_bounded_matrix log_rec_devs(1,n_ag,syr,nyr,-15.,15.,2);
	ofstream of4("iscam_rdev_mcmc.csv", ios::app);
	iter = 1;
	for(int yr = syr + sage; yr <= nyr; yr++){
	  if(iter == 1){
	    of4<<log_rec_devs(1)(yr);
	  }else{
	    of4<<","<<log_rec_devs(1)(yr);
	  }
	  iter++;
	}
	of4<<'\n';
	// output vulnerable biomass to all gears
	ofstream of5("iscam_vbt_mcmc.csv", ios::app);
	iter = 1;
	for(int fleet = 1; fleet <= nfleet; fleet++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(iter == 1){
	      of5<<vbt(1)(fleet)(yr);
	    }else{
	      of5<<","<<vbt(1)(fleet)(yr);
	    }
	    iter++;
	  }
	}
	of5<<'\n';
	// output fishing mortality as U (1-e^-F)
	ofstream of6("iscam_ut_mcmc.csv", ios::app);
	iter = 1;
	for(int fleet = 1; fleet <= nfleet; fleet++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(iter == 1){
	      of6<<1.0 - exp(-ft(1)(fleet)(yr));
	    }else{
	      of6<<","<<1.0 - exp(-ft(1)(fleet)(yr));
	    }
	    iter++;
	  }
	}
	of6<<'\n';
	// output natural mortality - (added for herring)
	// Assumes age-invariant M
	ofstream of7("iscam_m_mcmc.csv", ios::app);
	iter = 1;
	for(int sex = 1; sex <= n_gs; sex++){
	  for(int yr = syr; yr <= nyr; yr++){
	    if(iter == 1){
	      of7<<M(sex)(yr)(sage);
	    }else{
	      of7<<","<<M(sex)(yr)(sage);
	    }
	    iter++;
	  }
	}
	of7<<'\n';
	ofstream of8("iscam_index_fits_mcmc.csv", ios::app);
	iter = 1;
	for(int kk = 1; kk <= nItNobs; kk++){
	  for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	    if(iter == 1){
	      of8<<it_hat(kk,yr);
	    }else{
	      of8<<","<<it_hat(kk,yr);
	    }
	    iter++;
	  }
	}
	of8<<'\n';
	ofstream of9("iscam_index_residuals_mcmc.csv", ios::app);
	iter = 1;
	for(int kk = 1; kk <= nItNobs; kk++){
	  for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	    if(iter == 1){
	      of9<<epsilon(kk,yr);
	    }else{
	      of9<<","<<epsilon(kk,yr);
	    }
	    iter++;
	  }
	}
	of9<<'\n';
	ofstream of99("iscam_index_standardized_residuals_mcmc.csv", ios::app);
	iter = 1;
	for(int kk = 1; kk <= nItNobs; kk++){
	  for(int yr = 1; yr <= n_it_nobs(kk); yr++){
	    if(iter == 1){
	      of99<<std_epsilon(kk,yr);
	    }else{
	      of99<<","<<std_epsilon(kk,yr);
	    }
	    iter++;
	  }
	}
	of99<<'\n';
	ofstream of10("iscam_age_fits_mcmc.csv", ios::app);
        int gear, year, sex, a, row;
	for(k = 1; k <= nAgears; k++){
	  //of10<<"posterior"<<post_num<<"_gear"<<column(d3_A(k), -3)(1)<<"\n";
	  for(row = A_hat(k).rowmin(); row <= A_hat(k).rowmax(); row++){
	    gear = column(d3_A(k), -3)(1);
	    year = d3_A(k, row)(-5);
	    //gear = d3_A(k, row)(-3);
	    //area = d3_A(k, row)(-2);
	    //group = d3_A(k, row)(-1);
	    sex = d3_A(k, row)(0);
	    of10<<gear<<","<<post_num<<","<<year<<","<<sex<<",";
	    for(a = sage; a < nage; a++){
	      of10<<A_hat(k, row, age(a))<<",";
	    }
	    of10<<A_hat(k, row, age(nage))<<"\n";
	  }
	}
	ofstream of11("iscam_age_residuals_mcmc.csv", ios::app);
	for(k = 1; k <= nAgears; k++){
	  //of11<<"posterior"<<post_num<<"_gear"<<column(d3_A(k), -3)(1)<<"\n";
	  for(int row = A_nu(k).rowmin(); row <= A_nu(k).rowmax(); row++){
	    gear = column(d3_A(k), -3)(1);
	    year = d3_A(k, row)(-5);
	    //gear = d3_A(k, row)(-3);
	    //area = d3_A(k, row)(-2);
	    //group = d3_A(k, row)(-1);
	    sex = d3_A(k, row)(0);
	    of11<<gear<<","<<post_num<<","<<year<<","<<sex<<",";
	    for(a = sage; a < nage; a++){
	      of11<<A_nu(k, row, age(a))<<",";
	    }
	    of11<<A_nu(k, row, age(nage))<<"\n";
	  }
	}
	ofstream of12("iscam_selectivity_mcmc.csv", ios::app);
	int block;
	int last_block;
	for(k = 1; k <= ngear; k++){
	  last_block = n_sel_blocks(k);
	  for(block = 1; block <= last_block; block++){
	    if(block == 1){
	      if(block == last_block){
	        if(nsex == 2){
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<nyr<<
	            ",1,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<nyr<<
	            ",2,"<<exp(sel_par_m(k)(block)(1))<<","<<exp(sel_par_m(k)(block)(2))<<"\n";
		}else{
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<nyr<<
	            ",0,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
		}
	      }else{
	        if(nsex == 2){
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<sel_blocks(k, block+1)-1<<
	            ",1,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<sel_blocks(k, block+1)-1<<
	            ",2,"<<exp(sel_par_m(k)(block)(1))<<","<<exp(sel_par_m(k)(block)(2))<<"\n";
		}else{
	          of12<<k<<","<<post_num<<","<<block<<","<<syr<<","<<sel_blocks(k, block+1)-1<<
	            ",0,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
		}
	      }
	    }else{
	      if(block == last_block){
	        if(nsex == 2){
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<nyr<<
	            ",1,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<nyr<<
	            ",2,"<<exp(sel_par_m(k)(block)(1))<<","<<exp(sel_par_m(k)(block)(2))<<"\n";
		}else{
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<nyr<<
	            ",0,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
		}
	      }else{
	        if(nsex == 2){
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<sel_blocks(k, block+1)-1<<
	            ",1,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<sel_blocks(k, block+1)-1<<
	            ",2,"<<exp(sel_par_m(k)(block)(1))<<","<<exp(sel_par_m(k)(block)(2))<<"\n";
		}else{
	          of12<<k<<","<<post_num<<","<<block<<","<<sel_blocks(k, block)<<","<<sel_blocks(k, block+1)-1<<
	            ",0,"<<exp(sel_par_f(k)(block)(1))<<","<<exp(sel_par_f(k)(block)(2))<<"\n";
		}
	      }
	    }
	  }
	}

	ofs.flush();
	of1.flush();
	of2.flush();
	of3.flush();
	of4.flush();
	of5.flush();
	of6.flush();
	of7.flush();
	of8.flush();
	of9.flush();
	of10.flush();
	of11.flush();
	of12.flush();

	for(int ii = 1; ii <= n_tac; ii++){
	  LOG<<"TAC for projection #"<<ii<<": "<<tac(ii)<<"\n";
	  projection_model(tac(ii));
	}

	/*
	  This routine conducts population projections based on
	    the estimated values of theta.  Note that all variables
	    in this routine are data type variables.
	    Arguments:
	      tac is the total allowable catch that must be allocated
	      to each gear type based on dAllocation(k)
	      theta(1) = log_ro
	      theta(2) = h
	      theta(3) = log_m_female
	      if nsex == 1:
	      theta(4) = log_avgrec
	      theta(5) = log_recinit
	      theta(6) = rho
	      theta(7) = vartheta
	      if nsex == 2:
	      theta(4) = log_m_male
	      theta(5) = log_avgrec
	      theta(6) = log_recinit
	      theta(7) = rho
	      theta(8) = vartheta
	    NOTES:
	      - Projections are based on average natural mortality and fecundity.
	      - Selectivity is based on selectivity in terminal year.
	      - Average weight-at-age is based on mean weight in the last N years,
	      - where N is the difference in years given in items 3 and 4 of the projection
	      - file control options "years for average fecundity/weight-at-age in projections"
	*/
FUNCTION void projection_model(const double& tac);

	static int run_num = 0;
	run_num++;
	if(run_num == 1){
	  LOG<<"Running projection model...\n";
	}

	// n_projyr is defined in the .pfc file
	int pyr = nyr + n_projyr;
	int i, j, k, sex;
	double tau, phib, so, bo, beta, xx, rtt, et;
	BaranovCatchEquation cBaranov;

	dvector p_sbt(syr, pyr + 1);
	dvector p_rt(syr + sage, pyr);
	dmatrix lx(1, nsex, sage, nage);
	dmatrix lw(1, nsex, sage, nage);
	dmatrix m_bar(1, nsex, sage, nage);
	dmatrix fa_bar(1, nsex, sage, nage);
	dmatrix p_ct(1, nsex, 1, ngear);
	d3_array va_bar(1, nsex, 1, ngear, sage, nage);
	d3_array p_ft(1, nsex, nyr + 1, pyr + 1, 1, ngear);
	d3_array p_N(1, nsex, syr, pyr + 2, sage, nage);
	d3_array p_Z(1, nsex, syr, pyr + 1, sage, nage);

	p_sbt.initialize();
	p_rt.initialize();
	lx.initialize();
	lw.initialize();
	m_bar.initialize();
	fa_bar.initialize();
	p_ct.initialize();
	va_bar.initialize();
	p_ft.initialize();
	p_N.initialize();
	p_Z.initialize();

	// Fill arrays with historical values, no projections yet
	for(i = syr; i <= nyr; i++){
	  for(sex = 1; sex <= nsex; sex++){
	    p_N(sex, i) = value(N(sex, i));
	    // There is only one group, so sbt and rt have
	    // a hardwired 1 here
	    if(sex == 1){
	      p_sbt(i) = value(sbt(1, i));
	      if(i >= syr + sage){
	        p_rt(i) = value(rt(1, i));
	      }
	    }
	    p_Z(sex, i) = value(Z(sex)(i));
	  }
	}

	// Selectivity and gear allocation values
	for(k = 1; k <= ngear; k++){
	  for(sex = 1; sex <= nsex; sex++){
	    p_ct(sex, k) = dAllocation(k) * tac;
	    va_bar(sex, k) = exp(value(log_sel(k, sex, nyr)));
	  }
	}

 	tau = value(sqrt(1.0 - rho) * varphi);
	phib = 0;
	for(sex = 1; sex <= nsex; sex++){
	  // fa_bar is average fecundity by sex
	  fa_bar(sex) = elem_prod(dWt_bar(sex), ma(sex));
	  // ma_bar is average natural mortality by sex for years
	  // pf_cntrl(1) and pf_cntrl(2) which are in the .pfc file
	  m_bar(sex) = colsum(value(M(sex).sub(pf_cntrl(1), pf_cntrl(2))));
	  m_bar(sex) /= pf_cntrl(2) - pf_cntrl(1) + 1;

	  // Derive stock recruitment parameters
	  // lx is survivorship of spawning biomass by age
	  lx(sex, sage) = 1.0;
	  lw(sex, sage) = 1.0;
	  for(j = sage; j <= nage; j++){
	    if(j > sage){
	      lx(sex, j) = lx(sex, j - 1) * mfexp(-m_bar(sex, j - 1));
	    }
	    lw(sex, j) = lx(sex, j) * mfexp(-m_bar(sex, j) * d_iscamCntrl(13));
	  }
	  lw(sex, nage) /= 1.0 - mfexp(-m_bar(sex, nage));

	  phib += lw(sex) * fa_bar(sex);
	}
	so = value(kappa(1) / phib);
	bo = value(ro(1) * phib);
	beta = 1;
	switch(int(d_iscamCntrl(2))){
	  case 1: // Beverton-Holt
	    beta = value((kappa(1) - 1.) / bo);
	    break;
	  case 2: // Ricker
	    beta = value(log(kappa(1) / bo));
	    break;
	}

	// Simulate population into the future under constant tac policy
	for(i = nyr; i <= pyr + 1; i++){
	  for(sex = 1; sex <= nsex; sex++){
	    // ft(nyr) is a function of ct(nyr) not the tac so use ft(nyr)
	    // from the main model for nyr
	    if(i > nyr && i <= pyr + 1){
	      p_ft(sex, i) = cBaranov.getFishingMortality(
	        p_ct(sex),
	        m_bar(sex),
	        va_bar(sex),
	        p_N(sex, i),
	        dWt_bar(sex));
	      // Calculate total mortality in future years
	      p_Z(sex, i) = m_bar(sex);
	      for(k = 1; k <= ngear; k++){
	        p_Z(sex, i) += p_ft(sex, i, k) * va_bar(sex, k);
	      }
	    }
	    // Overwrite sbt(nyr) so that it does not include estimated Rt(nyr),
	    // which is highly uncertain. This will only be different from
	    // sbt(nyr) in the main model if recruits contribute to the spawning
	    // population, which is rare.
	    // d_iscamCntrl(13) is defined as: fraction of total mortality that
	    // takes place prior to spawning.
	    if(sex == 1){
	      // Overwrite the sbt the first time through the loop so the +=
	      // below doesn\'t add on to the estimated value
	      p_sbt(i) = 0;
	    }
	    p_sbt(i) +=
	      elem_prod(p_N(sex, i), mfexp(-p_Z(sex, i) * d_iscamCntrl(13))) * fa_bar(sex);
	    if(sex == 1){
	      // sage (age-1 typically) recruits with random deviate xx
	      // Note the random number seed is repeated for each tac level
	      // Note that this treatment of rec devs is different from
	      // historical model
	      xx = randn(nf + i) * tau;
	      if(i > nyr){
	        rtt = 1;
	        // Lagged spawning biomass
	        // (+1 because we want recruits for year i + 1)
	        et = p_sbt(i - sage + 1);
	        switch(int(d_iscamCntrl(2))){
	          case 1: // Beverton-Holt
	            rtt = (so * et / (1. + beta * et));
	            break;
 	          case 2: // Ricker
	            rtt = (so * et * exp(-beta * et));
	            break;
	        }
	        if(i <= pyr){
	          p_rt(i) = rtt;
	        }
	        // Next year\'s recruits
	        p_N(sex, i + 1, sage) = rtt * exp(xx - 0.5 * tau * tau);
	      }
	    }
	    // Update numbers at age in future years
	    // Next year\'s numbers
	    p_N(sex, i + 1)(sage + 1, nage) =
	      ++elem_prod(p_N(sex, i)(sage, nage - 1), exp(-p_Z(sex, i)(sage, nage - 1)));
	    p_N(sex, i + 1)(nage) +=
	      p_N(sex, i)(nage) * exp(-p_Z(sex, i, nage));
	  }
	}

	// write_proj_headers() and write_proj_output() are in include/utilities.h
	// and libs/utilities.cpp
	if(mceval_phase()){
	  ofstream ofsmcmc("iscam_mcmc_proj.csv", ios::app);
	  if(nf == 1 && run_num == 1){
	    write_proj_headers(ofsmcmc,
	                       syr,
	                       nyr,
	                       nfleet,
	                       n_ags,
	                       ngroup,
			       pyr,
	                       d_iscamCntrl(21),
			       d_iscamCntrl(13) && d_iscamCntrl(20));
	    ofsmcmc.flush();
	  }
	  write_proj_output(ofsmcmc,
	                    syr,
	                    nyr,
	                    nage,
	                    nfleet,
	                    n_ags,
	                    ngroup,
	                    tac,
	                    pyr,
	                    p_sbt,
	                    p_rt,
	                    p_ft,
	                    p_N,
	                    M,
	                    ma,
	                    dWt_bar,
	                    ft,
	                    value(sbo(1)),
	                    fmsy,
	                    bmsy,
	                    d_iscamCntrl(21),
	                    d_iscamCntrl(13) && d_iscamCntrl(20));
	  ofsmcmc.close();
	}else{
	  ofstream ofsmpd("iscam_mpd_proj.csv", ios::app);
	  //LOG<<"Running MPD projections, run_num = "<<run_num<<"\n";
	  if(run_num == 1){
	    write_proj_headers(ofsmpd,
	                       syr,
	                       nyr,
	                       nfleet,
	                       n_ags,
	                       ngroup,
			       pyr,
	                       !d_iscamCntrl(13),
	                       d_iscamCntrl(13) && d_iscamCntrl(20));
	    ofsmpd.flush();
	  }
	  //LOG<<"ft:"<<ft<<"\n";
	  //LOG<<"p_sbt:\n"<<p_sbt<<"\n";
	  //LOG<<"p_rt:\n"<<p_rt<<"\n";
	  //LOG<<"p_ft:\n"<<p_ft<<"\n\n";
	  write_proj_output(ofsmpd,
	                    syr,
	                    nyr,
	                    nage,
	                    nfleet,
	                    n_ags,
	                    ngroup,
	                    tac,
	                    pyr,
	                    p_sbt,
	                    p_rt,
	                    p_ft,
	                    p_N,
	                    M,
	                    ma,
	                    dWt_bar,
	                    ft,
	                    value(sbo(1)),
	                    fmsy,
	                    bmsy,
	                    !d_iscamCntrl(13),
	                    d_iscamCntrl(13) && d_iscamCntrl(20));
	  ofsmpd.flush();
	}
	if(!mceval_phase()){
	  //LOG<<"Finished projection model for TAC = "<<tac<<"\n";
	}

TOP_OF_MAIN_SECTION
	// These lines make all stdout and stderr go to the file
	int fd = open("admb_runtime.log", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
	dup2(fd, 1);
	dup2(fd, 2);
	close(fd);
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e8);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
	/*
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
	#include <unistd.h>
	#include <fcntl.h>
	#include "/usr/bin/admb/contrib/statslib/statsLib.h"
	#include "../../include/baranov.h"
	#include "../../include/ddirmultinom.h"
	#include "../../include/LogisticNormal.h"
	#include "../../include/LogisticStudentT.h"
	#include "../../include/msy.h"
	#include "../../include/slowmsy.h"
	#include "../../include/msy.hpp"
	#include "../../include/multinomial.h"
	#include "../../include/utilities.h"
	#include "../../include/Logger.h"
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;
	// Variables to store results from DIC calculations.
	double dicNoPar = 0;
	double dicValue = 0;

	/*
	  Slow msy routine to test ref points
	  Called by calcReferencePoints if turned on
	*/
FUNCTION void run_FRP()
	//Reference points
	if(last_phase()){
	  LOG<<"\n*********Getting reference points the slow way************\n";
	  LOG<<"*******************************************\n\n";
	}
	// ftest is a vector of F\'s from 0 to the number selected in the control
	// file with the number of elements necessary to match the precision
	// selected in the control file.
	int vec_size = (int)(d_iscamCntrl(19) / d_iscamCntrl(18)) + 1;
	dvector ftest(1, vec_size);
	ftest.fill_seqadd(0, d_iscamCntrl(18));
	int Nf = size_count(ftest);
	double Fmsy;
	double MSY;
	double Bmsy;
	double B0;
	// Matrix for putting numerically derived equilibrium catches for
	// Calculating MSY and FMSY (in R)
	dvector Ye(1,Nf);
	// Matrix for putting numerically derived equilibrium catches for
	// calculating MSY and FMSY (in R)
	dvector Be(1,Nf);
	dvector ye(1, Nf);
	dvector be(1, Nf);
	slow_msy(ftest, Ye, Be, MSY, Fmsy,
	         Bmsy, B0, sage, nage, nyr,
	         M, dWt_bar, ma, ro, kappa,
	         log_sel, d_iscamCntrl, pf_cntrl);
	if(d_iscamCntrl(13)){
	  // If the fraction of total mortality that takes place
	  //  before spawning is greater than 0, set the report
	  //  file outputs to the results from the slow_msy
	  //  routine
	  fmsy = Fmsy;
	  msy = MSY;
	  bmsy = Bmsy;
	  ye = Ye;
	  be = Be;
	  bo = B0;
	}
	ofstream ofsr("TEST_frp.rep");
	ofsr<<"Max F"<<'\n'<<d_iscamCntrl(19)<<'\n';
	ofsr<<"Precision for F"<<'\n'<<d_iscamCntrl(18)<<'\n';
	ofsr<<"Fmsy"<<'\n'<<Fmsy<<'\n';
	ofsr<<"MSY"<<'\n'<<MSY<<'\n';
	ofsr<<"Bmsy"<<'\n'<<Bmsy<<'\n';
	ofsr<<"ftest"<<'\n'<<ftest<<'\n';
	ofsr<<"Ye"<<'\n'<<Ye<<'\n';
	ofsr<<"Be"<<'\n'<<Be<<'\n';
	ofsr<<"B0"<<'\n'<<B0<<'\n';

FINAL_SECTION
	LOG<<"\n\nNumber of function evaluations: "<<nf<<'\n';
