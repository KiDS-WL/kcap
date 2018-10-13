///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - COSEBIs
*/
#include "COSEBIs.h"

//namespace COSEBIS_ {
extern "C" {
	const string sectionName = "cosebis";
	static double staticThing=1;
	//const char *paramsSection = "cosebis_params";
	const string shear_cl = SHEAR_CL_SECTION;
	const number MinPowerCOSEBIs=0.1;
	const number MaxPowerCOSEBIs=1e6;
	const int PowerTableNumberCOSEBIs=200;

	///define a structure with everything that is needed to be read in setup and sent to execute
	// Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 

	struct COSEBIs_config {
		string input_section_name;
		string output_section_name;
		int n_max;
		number theta_min;
		number theta_max;
		int input_option;//this can be either 0 for Cl or 1 for xi_pm
		int noisy_input;//This is either 0: no noise or 1: noisy
		COSEBIs *cosebis;
	} ;
	COSEBIs_config config;

	

  	void *setup(cosmosis::DataBlock *options) 
  	{
  		//options reads the ini file
  		//define config here and then read from options the relevant input quantities
  		clog<<"?????????????????????????????????????????????"<<endl;
  		clog<<"I'm in status in COSEBIs_cosmosis_inteface"<<endl;
  		clog<<"????????????????????????????????????????????"<<endl;
  
  		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    	const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
    	
  		status=options->get_val<string>(sectionName, string("output_section_name"), string("COSEBIs_results"), config.output_section_name);
  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config.theta_min);
  		if (status) 
		{
			clog<<"Could not load theta_min to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<1.<<endl;
			//exit(1);
		}
		else
		{
			clog<<"got the value of theta_min="<<config.theta_min<<endl;
		}
  		//if (status) return failure;
  		//clog<<"reading the value of theta_max"<<endl;
	    status=options->get_val<number>(sectionName, string("theta_max"),100., config.theta_max);
	    if (status) 
		{
			clog<<"Could not load theta_max to COSEBIs"<<endl;
			//exit(1);
		}
		else
		{
			clog<<"got the value of theta_max="<<config.theta_max<<endl;
		}
	    status=options->get_val<int>(sectionName, string("n_max"),10, config.n_max);
	    if (status) 
		{
			clog<<"Could not load n_max to COSEBIs"<<endl;
			//exit(1);
		}
		else
		{
			clog<<"got the value of nax_max"<<endl;
		}
		//This is either 0: Cl as input or 1: xi_+- as input
		//default is Cl. 
		status=options->get_val<int>(sectionName, string("input_option"), 0 , config.input_option);
		status=options->get_val<int>(sectionName, string("noisy_input"), 0 , config.noisy_input);
		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
	    cosebis->initialize(config.n_max,config.theta_min,config.theta_max,1);//npair set to one for now, will be set seperately in execute to the correct value
	    clog<<"checking which input option is chosen the default is Cl"<<endl;
		if(config.input_option==0)
		{
			clog<<"input is set to Cl"<<endl;
			cosebis->setWns(config.n_max);
			status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_cl"), config.input_section_name);
			if (status) 
				clog<<"Could not load n_max to COSEBIs"<<endl;
			else
				clog<<"got the value of nax_max"<<endl;
		}
		else if (config.input_option==1)
		{
			clog<<"input is set to 2PCFs"<<endl;
			cosebis->setTs(config.n_max);
			status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_shear"), config.input_section_name);
			if (status) 
				clog<<"Could not load n_max to COSEBIs"<<endl;
			else
				clog<<"got the value of nax_max"<<endl;
		}
		else
		{
			clog<<"ERROR: not a known input for COSEBIs, options are 0:Cl 1:xi_pm"<<endl;
			exit(1);
		}
	    
	    config.cosebis=cosebis;
  		return &config;
  		// config is sent to execute 
	}

	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
	    // Config is whatever you returned from setup above
	    // Block is the collection of parameters and calculations for
	    // this set of cosmological parameters
	    
	    COSEBIs_config *config= (COSEBIs_config*) config_in;
	    
	    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
	    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

	    //cout<<"config->"
	    clog<<"/////////////////\\\\\\\\\\\\\\\\\\\\\\"<<endl;
	    clog<<"I'm in execute in COSEBIs_cosmosis_inteface"<<endl;
	    clog<<"/////////////////\\\\\\\\\\\\\\\\\\\\\\"<<endl;
	    //clog<<"config->input_section_name="<<(config->input_section_name)<<endl;
	    //cout <<"type of config is: "<< typeid(config).name() << '\n';
	    //check if there is more than one input number of bins: nbin_a, nbin_b
	    
	    int num_z_bin_A;
		int num_z_bin_B;
	    //if (cosmosis::has_val(config->input_section_name, string("nbin_a")))
	    //{
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
		else
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
		}
		//clog<<"config->input_section_name="<<(config->input_section_name)<<endl;
		//clog<<"num_z_bin_A="<<num_z_bin_A<<endl;
		//clog<<"num_z_bin_B="<<num_z_bin_B<<endl;
		
		
		matrix En_mat;
		//does it have to be this format?
		if(config->input_option==0)
		{
			vector<number> ell,logell;
			status = block->get_val(config->input_section_name, string("ell"), ell);
			int nell=ell.size();
			for(int i=0; i<nell; i++)
			{
				logell.push_back(log(ell[i]));
			}
			//clog<<ell[0]<<'\t'<<ell[1]<<'\t'<<ell[nell-1]<<endl;
			if ((ell[0]>MinPowerCOSEBIs) || (ell[nell-1]<MaxPowerCOSEBIs) || (nell<PowerTableNumberCOSEBIs))
			{
				clog<<"*************************************"<<endl;
				clog<<"****************WARNING**************"<<endl;
				clog<<"ell range or number of points not sufficient for COSEBIs"<<endl;
				clog<<"your input ell_min="<<ell[0]<<"  ell_max="<<ell[nell-1]<<" n_ell="<<nell<<endl;
				clog<<" For a higher precision at least use this range:"<<endl;
				clog<<"ell_min="<<MinPowerCOSEBIs<<"  ell_max="
					<<MaxPowerCOSEBIs<<" n_ell="<<PowerTableNumberCOSEBIs<<endl;
				clog<<"*************END OF WARNING***********"<<endl;
				clog<<"*************************************"<<endl;
			}
			
			if (status) 
			{
				clog<<"Could not load ell in C_ell to COSEBIs"<<endl;
				return status;
			}

			vector <vector<number> > InputPower_vec_vec;
			int nPairs=0;
			for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
				{
					// read in C(l)
					vector<number> C_ell;
					string name_in="bin_"+toString(j_bin)+"_"+toString(i_bin);
					clog<<"loading Cl into InputPower_vec_vec, input name is: "<<name_in<<endl;
		    		status = block->get_val(config->input_section_name, name_in, C_ell);
					if (status) 
					{
						clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in C_ell to COSEBIs"<<endl;
						return status;
					}
					
					InputPower_vec_vec.push_back(C_ell);
					///put Cl in a vector to send to COSEBIs.
					nPairs++;
				}
			}
			config->cosebis->setZbins(nPairs);
			config->cosebis->setPower(logell,InputPower_vec_vec);
			En_mat=config->cosebis->calEn();
			int n_max=config->n_max;
		}
		///fix
		else if (config->input_option==1)
		{
			vector<number> theta, logTheta;
			status = block->get_val(config->input_section_name, string("theta"), theta);
			if (status) 
			{
				clog<<"Could not load theta in 2PCFs to COSEBIs"<<endl;
				return status;
			}
			int nTheta=theta.size();
			for(int i=0; i<nTheta; i++)
				logTheta.push_back(log(theta[i]));

			vector <vector<number> > xiP_vec_vec;
			vector <vector<number> > xiM_vec_vec;
			int nPairs=0;
			for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
				{
					vector<number> xiP_vec,xiM_vec;
					string name_in="xiplus_"+toString(j_bin)+"_"+toString(i_bin);
		    		status = block->get_val(config->input_section_name, name_in, xiP_vec);
					if (status) 
					{
						clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in XiP to COSEBIs"<<endl;
						return status;
					}
					string name_in="ximinus_"+toString(j_bin)+"_"+toString(i_bin);
		    		status = block->get_val(config->input_section_name, name_in, xiM_vec);
					if (status) 
					{
						clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in XiM to COSEBIs"<<endl;
						return status;
					}
					xiP_vec_vec.push_back(xiP_vec);
					xiM_vec_vec.push_back(xiM_vec);
					nPairs++;
				}
			}
			config->cosebis->setZbins(nPairs);
			config->cosebis->setKsi(logTheta,xiP_vec_vec,xiM_vec_vec);
			if(config->noisy_input==0)
				En_mat=config->cosebis->calEn();
			else if (config->noisy_input==1)
				En_mat=config->cosebis->calEn();
			else
			{
				clog<<"noisy_input not a recognized option, use either 0 for no noise or 1 for noisy"<<endl;
				exit(1);
			}
		}
		else
		{
			clog<<"input_option not recognized use either 0 for Cl or 1 for 2PCFs"<<endl;
			exit(1);
		}

		int n_max=config->n_max;
		vector<number> En_vec(n_max);
		vector<int> n_vals(n_max);
		int p1=0;
		string name_En;
		for(int i_bin=0; i_bin<num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<num_z_bin_B; j_bin++) 
			{
				int m=0;
				//int p1=cosebis.calP(nBins,bin1,bin2);
				for(int n1=n_max*p1,m=0 ;n1<n_max*(p1+1) ;n1++,m++)
					En_vec[m]=En_mat.get(n1);
				if(config->input_option==0)
					name_En=string("cosebis_EnFromCl_bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
				else if (config->input_option==1)
					name_En=string("cosebis_EnFrom2PCFs_bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
				status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
				p1++;
			}
		}

		for(int n=0;n<n_max;n++)
			n_vals[n]=n+1;

	    status = block->put_val<vector<int> >(config->output_section_name, string("cosebis_n"), n_vals);

	    ///Calculate likelihood here
	    number likelihood_val=2.+staticThing; //just testing now
	    staticThing+=0.01;
	    clog<<"likelihood="<<likelihood_val<<endl;
	    status = block->put_val<number>(LIKELIHOODS_SECTION, "COSEBIs_like", likelihood_val);
	    return status;
	}
}// end of extern C

//TODO: 
/*
1- make three cosebis interfaces for one point measurement
1.1- from Cl
1.2 from Ksi
1.3 from noisy Ksi
2- make one for likelihood analysis
2.1- measure the covariance matrix, if it is already calculated dont calculate again unless a new one is needed for the new sets of 
parameters
2.2- put an option to read from a file
2.3- another option to calculate it from input xi_+- files (sims)
3- test the likellihood



    
