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
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 

	struct COSEBIs_config {
		string input_section_name;
		string output_section_name;
		int n_max;
		number theta_min;
		number theta_max;
		int IsItBmodes;
		COSEBIs *cosebis;
	} ;
	COSEBIs_config config;

	

  	void *setup(cosmosis::DataBlock *options) 
  	{
  		//options reads the ini file
  		//define config here and then read from options the relevant input quantities
  		clog<<"?????????????????????????????????????????????"<<endl;
  		clog<<"I'm in status in COSEBIs_Cl_cosmosis_inteface"<<endl;
  		clog<<"????????????????????????????????????????????"<<endl;
  
  		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    	const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
    	
  		status=options->get_val<string>(sectionName, string("output_section_name"), string("COSEBIs_results"), config.output_section_name);
  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config.theta_min);
  		if (status) 
		{
			clog<<"Could not load theta_min to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<1.<<endl;
		}
		else
			clog<<"got the value of theta_min="<<config.theta_min<<endl;
	    status=options->get_val<number>(sectionName, string("theta_max"),100., config.theta_max);

	    if (status) 
			clog<<"Could not load theta_max to COSEBIs"<<endl;
		else
			clog<<"got the value of theta_max="<<config.theta_max<<endl;

	    status=options->get_val<int>(sectionName, string("n_max"),10, config.n_max);
	    if (status) 
			clog<<"Could not load n_max to COSEBIs"<<endl;
		else
			clog<<"got the value of n_max="<<config.n_max<<endl;
		//   initialize COSEBIs
		status=options->get_val<int>(sectionName, string("is_it_bmodes"),0, config.IsItBmodes);
		if(config.IsItBmodes)
		{
			clog<<"We are going to calculate B-modes"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!NOTE: make sure that the input Cl is for B-modes !!!!!!!"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		}
		else
			clog<<"Caluclating E-modes"<<endl;

		string WnFolderName,TnFolderName,OutputTnFolderName;
		bool defaultFolderName=false;
		status=options->get_val<string>(sectionName, string("Wn_Output_FolderName"), WnFolderName);
		if(status)
		{
			clog<<"Could not find WnLog folder name in Wn_Output_FolderName, setting to default"<<endl;
			defaultFolderName=true;
		}
		else
			clog<<"WnLog folder name is:"<<WnFolderName<<endl;

		status=options->get_val<string>(sectionName, string("Roots_n_Norms_FolderName"), TnFolderName);
		if(status)
		{
			clog<<"Could not find Root and Norms folder name in Roots_n_Norms_FolderName, setting to default"<<endl;
			defaultFolderName=true;
		}
		else
			clog<<"Root and Norms folder name is:"<<TnFolderName<<endl;


		status=options->get_val<string>(sectionName, string("Tn_Output_FolderName"), OutputTnFolderName);
		if(status)
		{
			clog<<"Could not find T_pm folder name in Tn_Output_FolderName, setting to default"<<endl;
			defaultFolderName=true;
	  	}
		else
			clog<<"T_pm folder name is:"<<OutputTnFolderName<<endl;

		COSEBIs *cosebis = new COSEBIs();
		
		if(defaultFolderName)
			cosebis->initialize(config.n_max,config.theta_min,config.theta_max,1);//npair set to one for now, will be set seperately in execute to the correct value
		else
			cosebis->initialize(config.n_max,config.theta_min,config.theta_max,1 //npair set to one for now, will be set seperately in execute to the correct value
				,TnFolderName,WnFolderName,OutputTnFolderName);

		cosebis->setWns(config.n_max);
		status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_cl"),
			config.input_section_name);
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
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
		else
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
		}
		matrix En_mat;
		vector<number> ell,logell;
		status = block->get_val(config->input_section_name, string("ell"), ell);
		int nell=ell.size();
		for(int i=0; i<nell; i++)
			logell.push_back(log(ell[i]));
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
				//clog<<name_in<<endl;
				status = block->get_val(config->input_section_name, name_in, C_ell);
				string name_Cl=string("Cl_bin_")+toString(j_bin)+string("_")+toString(i_bin)+string(".ascii");
				//status = block->put_val<vector<number> >(config->output_section_name, name_Cl, C_ell);
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
		//status = block->put_val<vector<number> >(config->output_section_name, string("ell.ascii"), ell);
		config->cosebis->setZbins(nPairs);
		config->cosebis->setPower(logell,InputPower_vec_vec);
		clog<<"ell="<<15<<"    Power="<<config->cosebis->ReturnPower(15.,0);
		En_mat=config->cosebis->calEn();
		// En_mat.printOut("")
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
				name_En=string("bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
				status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
				p1++;
			}
		}
		status = block->put_val<int>(config->output_section_name, string("n_mode"), n_max);
		status = block->put_val<bool>(config->output_section_name, string("b_modes"), config->IsItBmodes);
		status = block->put_val<double>(config->output_section_name, string("theta_min"), config->theta_min);
		status = block->put_val<double>(config->output_section_name, string("theta_max"), config->theta_max);

		for(int n=0;n<n_max;n++)
			n_vals[n]=n+1;

	    status = block->put_val<vector<int> >(config->output_section_name, string("cosebis_n"), n_vals);
	    //string EnfileName="EnTH_Cosmosis_4bins_KiDSBestfit_0.50-100.00";
	    //En_mat.store((EnfileName+string(".bin")).c_str());
	    //En_mat.printOut((EnfileName+string(".ascii")).c_str(),10);
	    return status;
	}
}// end of extern C

//TODO: 
/*
1.2 from Ksi
1.3 from noisy Ksi
2.1- measure the covariance matrix, if it is already calculated dont calculate again unless a new one is needed for the new sets of 
parameters
2.2- put an option to read from a file
2.3- another option to calculate it from input xi_+- files (sims)
3- test the likellihood
*/


    
