///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - BandPower
*/
#include "BandPower_W.h"


extern "C" 
{
	const string shear_cl = SHEAR_CL_SECTION;
	const int NLBINS=5000;
	const number LLOW=0.1;
	const number LHIGH= 1e4;


	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as BandPower. BandPower is a class that produces En/Bn and 
	//their covariance, etc. 

	int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
	{
	
		auto status = options->get_val(OPTION_SECTION, name, parameter);
		if (status!=DBS_SUCCESS) 
		{
			parameter = "";
			cerr<< "Could not find or understand parameter in BandPower section: " 
				<< name << std::endl; 
			return 1;
		}
		return 0;
	}

	struct BandPower_config {
		string sectionName;
		string input_section_name;
		string output_section_name;
		number theta_min;
		number theta_max;
		int IsItBmodes;
		bool calCov;
		bool calNoiseCov;
		string Cov_En_name;
		int nBins; //this is only needed if the noise only covariance is to be estimated from input nPair
		BandPower_W *bandpower;
		number sigma_m;
	} ;
	BandPower_config config;

	

  	void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
  	{
  		config.sectionName=OPTION_SECTION;
  		string sectionName=config.sectionName;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
		
		// This is where the outputs will be saved
		status=options->get_val<string>(sectionName, string("output_section_name"), config.output_section_name);
		if (status) 
		{
			clog<<"Could not load out_section_name to BandPower, ";
			clog<<"setting to default: band_power"<<endl;
			config.output_section_name=string("band_power");
		}
		else
			clog<<"got the value of output_section_name:"<<config.output_section_name<<endl;

		//get input section name, default= shear_cl
		status=options->get_val<string>(sectionName, string("input_section_name"), config.input_section_name);
		if (status) 
		{
			clog<<"Could not load input_section_name to band power,";
			clog<<" setting to default: shear_cl"<<endl;
			config.input_section_name=string("shear_cl");
		}
		else
			clog<<"got the value of input_section_name:"<<config.input_section_name<<endl;
	
		//minimum theta value that is used to estimate band powers from xi_pm, in arcmins
  		status=options->get_val<number>(sectionName, string("theta_min") , config.theta_min);
  		if (status) 
		{
			clog<<"Could not load theta_min to BandPower"<<endl;
			clog<<"setting it to the default value of 1 arcmins"<<endl;
			config.theta_min=1.;
		}
		else
			clog<<"got the value of theta_min="<<config.theta_min<<endl;

		//maximum theta value that is used to estimate band powers from xi_pm, in arcmins
	    status=options->get_val<number>(sectionName, string("theta_max"), config.theta_max);
	    if (status) 
	    {
			clog<<"Could not load theta_max to BandPower, setting it to the default value of 100 arcmin"<<endl;
			config.theta_max=100.;
	    }
		else
			clog<<"got the value of theta_max="<<config.theta_max<<endl;

		//first lets look for a file with l_min l_max for each band power bin
		string l_min_max_file;
		//vectors that contain the min and max value of ell for each band power bin
		vector<number> l_min_vec,l_max_vec;
		status=options->get_val<string>(sectionName, string("l_min_max_file"), l_min_max_file);
		int nBands;
		if (status)
		{
			number l_min,l_max;
			clog<<"Could not load l_min_max_file to BandPower"; 
			clog<<"going to look for nBands and l_min and l_max instead"<<endl;
			//number of bands
		    status=options->get_val<int>(sectionName, string("nBands"), nBands);
		    if (status) 
				clog<<"Could not load n_max to BandPower, setting it to the default value of 10"<<endl;
			else
				clog<<"got the value of nBands="<<nBands<<endl;

			status=options->get_val<number>(sectionName, string("l_min"), l_min);
			if (status) 
		    {
				clog<<"Could not load l_min to BandPower, setting it to 100."<<endl;
				l_min=100.;
		    }
			else
				clog<<"got the value of l_min="<<l_min<<endl;

			status=options->get_val<number>(sectionName, string("l_max"), l_max);
			if (status) 
		    {
				clog<<"Could not load l_max to BandPower, setting it to 1000"<<endl;
				l_max=1000.;
		    }
			else
				clog<<"got the value of l_max="<<l_max<<endl;
			
			for(int b=0; b<nBands; b++)
			{
				number lmin=exp(log(l_min)+log(l_max/l_min)/(nBands)*b);
				number lmax=exp(log(l_min)+log(l_max/l_min)/(nBands)*(b+1));
				l_min_vec.push_back(lmin);
				l_max_vec.push_back(lmax);
				clog<<"b="<<b<<" lmin="<<lmin<<" lmax="<<lmax<<endl;
			}
		}
		else
		{
			matrix l_min_max_mat;
			l_min_max_mat.readFromASCII_marika(l_min_max_file.c_str());
			nBands=l_min_max_mat.rows;
			for(int b=0; b<nBands; b++)
			{
				l_min_vec.push_back(l_min_max_mat.get(0,b));
				l_max_vec.push_back(l_min_max_mat.get(1,b));
			}
		}

		//If set to 1 then calculates B-mode band powers, needs b-mode Cls as input.s
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


		string FolderName;
		bool defaultFolderName=false;
		status=options->get_val<string>(sectionName, string("Output_FolderName"), FolderName);
		if(status)
		{
			clog<<"Could not find the folder name for weight functions in W_Output_FolderName, setting to default"<<endl;
			defaultFolderName=true;
		}
		else
			clog<<"W folder name is:"<<FolderName<<endl;


		number  Delta_x;
		status=options->get_val<number>(sectionName, string("Delta_x"), Delta_x);
		if(status)
		{
			clog<<"Could not find Delta_x, setting to default:"<<endl;
			Delta_x=0.5; //Need to check this
		}
		else
			clog<<"Delta_theta is:"<<Delta_x<<endl;

		//number Delta_x=log(Delta_theta*arcmin);

		string Response_function_type;
		status=options->get_val<string>(sectionName, string("Response_function_type"), Response_function_type);
		if(status)
		{
			clog<<"Could not find the response function type setting to default: tophat"<<endl;
			Response_function_type="tophat";
		}
		else
			clog<<"response function type is:"<<Response_function_type<<endl;

		bool Analytic=true; //set this to be read from the ini file
		int Analytic_int;
		status=options->get_val<int>(sectionName, string("Analytic"), Analytic_int);
		if(status)
		{
			clog<<"Could not find Analytic, setting to default: true"<<endl;
			//Analytic=0.5; //Need to check this
		}
		else
		{
			if(Analytic_int==0)
			{
				Analytic=false;
				clog<<"setting analytic to false"<<endl;
			}
			else
				clog<<"setting analytic to true"<<endl;
		}
		//   initialize BandPower
		//BandPower_W *W = new BandPower_W();
		// need to think about this one
		// int bessel_order=0;
		// bool noApodise=true;
		// clog<<"going to initialize W "<<endl;
		// //now testing W
		// //Delta_x=0.5;
		// W->initialize(config.theta_min*arcmin,config.theta_max*arcmin,Response_function_type,
		// 		l_min_vec,l_max_vec,bessel_order,noApodise,Delta_x,Analytic,LLOW,LHIGH,NLBINS);
		// FolderName="../cosmosis_in_out/BandPower_results/";
		// W->setBandPower_WName(FolderName,"g","W_noAp");
		// for(int b=0; b<nBands; b++)
		// {
		// 	//W->print_integrant(noApodise,0,1000,b);
		// 	W->set(b,0);
		// 	W->set(b,2);
		// 	W->set(b,4);
		// }
		// W->setBandPower_WName(FolderName,"g","W_Ap");
		// noApodise=false;
		// W->set_noApodise(noApodise);
		// //clog<<"!!!!!!!!!!!!!!!!!!!setting Ap_more!!!!!!!!!!!!!!!!!!"<<endl;
		// for(int b=0; b<nBands; b++)
		// {
		// 	//W->print_integrant(noApodise,0,1000,b);
		// 	W->set(b,0);
		// 	//exit(1);
		// 	W->set(b,2);
		// 	W->set(b,4);
		// }
		BandPower *BP0 = new BandPower();
		BandPower *BP2 = new BandPower();
		BandPower *BP4 = new BandPower();
		BP0->initialize();
		
  		return &config;
  		// config is sent to execute 
	}

	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
		// Config is whatever you returned from setup above
		// Block is the collection of parameters and calculations for
		// this set of cosmological parameters
		
		BandPower_config *config= (BandPower_config*) config_in;
		
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//get cl from cosmosis

		//calculate BP

		//save BP

		

	    return status;
	}
}// end of extern C


    
