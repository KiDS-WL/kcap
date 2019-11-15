///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - COSEBIs
*/
#include "COSEBIs.h"
#include "calChiS.h"


extern "C" {
	//const string sectionName = "cosebis";
	const string shear_cl = SHEAR_CL_SECTION;
	const number MinPowerCOSEBIs=0.1;
	const number MaxPowerCOSEBIs=1e6;
	const int PowerTableNumberCOSEBIs=200;
	const string cterm_section= "shear_c_bias";

	typedef struct COSEBIs_config {
		string sectionName;
		string input_section_name;
		string output_section_name;
		int n_max;
		int IsItBmodes;
		number theta_min;
		number theta_max;
		int input_option;//this can be either 0 for Cl or 1 for xi_pm
		COSEBIs *cosebis;
		matrix En_data;
		matrix Cov_mat;
		bool add_cterm;
		matrix En_cos4phi;
		matrix En_sin4phi;
		bool add_2D_cterm;
		matrix En_2D;
		bool dolikelihood;
		//vector<number> m_bias;
		//bool apply_m_bias;
	} COSEBIs_config;


	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 
	int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
	{
	
		auto status = options->get_val(OPTION_SECTION, name, parameter);
		if (status!=DBS_SUCCESS) 
		{
			parameter = "";
			cerr<< "Could not find or understand parameter in cosebis section: " 
				<< name << std::endl; 
			return 1;
		}
		return 0;
	}
	
	void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
	//void *setup(cosmosis::DataBlock *options) 
	{
		//options reads the ini file
		//define config here and then read from options the relevant input quantities
		COSEBIs_config * config = new COSEBIs_config;

  		string sectionName=OPTION_SECTION;
  		config->sectionName=sectionName;

		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		config->add_cterm=false;
	
		//get output_section_name, default= "COSEBIs_results"
		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);
		if (status) 
		{
			clog<<"Could not load out_section_name to COSEBIs, ";
			clog<<"setting to default: cosebis_results"<<endl;
			config->output_section_name=string("cosebis");
		}
		else
			clog<<"got the value of output_section_name:"<<config->output_section_name<<endl;

		//get input section name, default= shear_cl
		status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
		if (status) 
		{
			clog<<"Could not load input_section_name to COSEBIs,";
			clog<<"setting to default: shear_cl"<<endl;
			config->input_section_name=string("shear_cl");
		}
		else
			clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;

		status=options->get_val<int>(sectionName, string("is_it_bmodes"),0, config->IsItBmodes);
		if(config->IsItBmodes)
		{
			clog<<"We are going to calculate B-modes"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!NOTE: make sure that the input Cl is for B-modes !!!!!!!"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		}
		else
			clog<<"Caluclating E-modes"<<endl;

		//get theta_min value
		status=options->get_val<number>(sectionName, string("theta_min"), config->theta_min);
		if (status) 
		{
			clog<<"Could not load theta_min to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<1.<<endl;
			config->theta_min=1.;
		}
		else
			clog<<"got the value of theta_min="<<config->theta_min<<endl;

		// status=options->get_val<vector<number>>(sectionName, string("m_bias_value"), config->m_bias);
		// if (status) 
		// {
		// 	clog<<"Could not load m_bias value to COSEBIs"<<endl;
		// 	clog<<"Not going to apply any corrections to the data vector"<<endl;
		// 	config->apply_m_bias=false;	
		// }
		// else
		// {
		// 	for(int i=0; i< config->m_bias.size(); i++)
		// 		clog<<"got the value of m_bias="<<config->m_bias[i]<<endl;
		// 	config->apply_m_bias=true;
		// }


		//get theta_max value
		status=options->get_val<number>(sectionName, string("theta_max"), config->theta_max);
		if (status) 
		{
			clog<<"Could not load theta_max to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<100.<<endl;
			config->theta_max=100.;
		}
		else
			clog<<"got the value of theta_max="<<config->theta_max<<endl;

		//get n_max value
		status=options->get_val<int>(sectionName, string("n_max"), config->n_max);
		if (status) 
		{
			clog<<"Could not load n_max to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<5<<endl;
			config->n_max=10;
		}
		else
			clog<<"got the value of n_max="<<config->n_max<<endl;

		//get input COSEBIs file name
		string InputCOSEBIsFileName;
		status=options->get_val<string>(sectionName, string("input_cosebis_filename"), InputCOSEBIsFileName);
		if (status) 
		{
			clog<<"Could not find the name of the input COSEBIs file, not going to calculate likelihood"<<endl;
			config->dolikelihood=false;
		}
		else
		{
			clog<<"got the name of the input COSEBIs file="<<InputCOSEBIsFileName<<endl;
			config->dolikelihood=true;
			///read En_data from file
			config->En_data.readFromASCII_marika((InputCOSEBIsFileName).c_str());
			//get input covariance file name
			string InputCovarianceFileName;
			status=options->get_val<string>(sectionName, string("input_covariance_filename"), InputCovarianceFileName);
			if (status) 
			{
				clog<<"Could not find the name of the input Covariance for COSEBIs file, no going to calculate likelihood"<<endl;
				config->dolikelihood=false;
			}
			else
			{
				clog<<"got the name of the input Covariance for COSEBIs file="<<InputCovarianceFileName<<endl;
				config->dolikelihood=true;
				///read Covariance from file
        		config->Cov_mat.readFromASCII_marika((InputCovarianceFileName).c_str());
        				string InputNGCovarianceFileName;
				status=options->get_val<string>(sectionName, string("input_NonGaussian_covariance_filename"), InputNGCovarianceFileName);
				if (status) 
					clog<<"Could not find the name of the input Non Gaussian Covariance for COSEBIs file"<<endl;
				else
				{
					clog<<"got the name of the input Non Gaussian Covariance for COSEBIs file="<<InputNGCovarianceFileName<<endl;
					matrix NGCov;
					NGCov.readFromASCII_marika((InputNGCovarianceFileName).c_str());
					config->Cov_mat+=NGCov;
				}


				string InputExtraCovarianceFileName;
				status=options->get_val<string>(sectionName, string("input_Extra_covariance_filename"), InputExtraCovarianceFileName);
				if (status) 
					clog<<"Could not find the name of the input Extra Covariance for COSEBIs file"<<endl;
				else
				{
					clog<<"got the name of the input Extra Covariance for COSEBIs file="<<InputExtraCovarianceFileName<<endl;
					clog<<"going to add it to the other covariances given"<<endl;
					matrix Cov;
					Cov.readFromASCII_marika((InputExtraCovarianceFileName).c_str());
					config->Cov_mat+=Cov;
				}
			}
		}

	
		string input_cterm_cos4phi_filename,input_cterm_sin4phi_filename;
		//get input cos4phi cterm file name
		status=options->get_val<string>(sectionName, string("input_cterm_cos4phi_filename"), input_cterm_cos4phi_filename);
		if (status) 
		{
			clog<<"Could not find the name of the input_cterm_cos4phi_filename for COSEBIs file"<<endl;
			config->add_cterm=false;
		}
		else
		{
			clog<<"got the name of the input_cterm_cos4phi_filename for COSEBIs file="<<input_cterm_cos4phi_filename<<endl;
			///read En_cos4phi from file
			config->En_cos4phi.readFromASCII_marika((input_cterm_cos4phi_filename).c_str());
			config->add_cterm=true;
		}

		//get input cos4phi cterm file name
		status=options->get_val<string>(sectionName, string("input_cterm_sin4phi_filename"), input_cterm_sin4phi_filename);
		if (status)
		{ 
			clog<<"Could not find the name of the input_cterm_sin4phi_filename for COSEBIs file"<<endl;
			config->add_cterm=false;
		}
		else
		{
			clog<<"got the name of the input_cterm_sin4phi_filename for COSEBIs file="<<input_cterm_sin4phi_filename<<endl;
			///read En_sin4phi from file
			config->En_sin4phi.readFromASCII_marika((input_cterm_sin4phi_filename).c_str());
		}

		if(config->add_cterm)
		{
			clog<<"going to add the constant c-term effect in the analysis"<<endl;
		}

		//get input 2D cterm file name
		string input_2Dcterm_filename;
		status=options->get_val<string>(sectionName, string("input_2Dcterm_filename"), input_2Dcterm_filename);
		if (status)
		{ 
			clog<<"Could not find the name of the input_2Dcterm_filename for COSEBIs file"<<endl;
			config->add_2D_cterm=false;
		}
		else
		{
			clog<<"got the name of the input_2Dcterm_filename for COSEBIs file="<<input_2Dcterm_filename<<endl;
			///read En_sin4phi from file
			config->En_2D.readFromASCII_marika((input_2Dcterm_filename).c_str());
			config->add_2D_cterm=true;
		}



		//get Wn, Tn and output Tn folder names
		string WnFolderName,TnFolderName,OutputTnFolderName;
		bool defaultFolderName=false;
		status=options->get_val<string>(sectionName, string("Wn_Output_FolderName"), WnFolderName);
		if(status)
		{
			clog<<"Could not find WnLog folder name in Wn_Output_FolderName, ";
			clog<<"setting to default: ./cosebis/WnLog/"<<endl;
			defaultFolderName=true;
		}
		else
			clog<<"WnLog folder name is:"<<WnFolderName<<endl;

		status=options->get_val<string>(sectionName, string("Roots_n_Norms_FolderName"), TnFolderName);
		if(status)
		{
			clog<<"Could not find Root and Norms folder name in Roots_n_Norms_FolderName,"; 
			clog<<" setting to default: ./cosebis/TLogsRootsAndNorms"<<endl;
			defaultFolderName=true;
		}
		else
			clog<<"Root and Norms folder name is:"<<TnFolderName<<endl;


		status=options->get_val<string>(sectionName, string("Tn_Output_FolderName"), OutputTnFolderName);
		if(status)
		{
			clog<<"Could not find T_pm folder name in Tn_Output_FolderName, ";
			clog<<"setting to default: ./cosebis/TpnLog/"<<endl;
			defaultFolderName=true;
	  	}
		else
			clog<<"T_pm folder name is:"<<OutputTnFolderName<<endl;

		


		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
		if(defaultFolderName)
		{
			clog<<"going to initialize COSEBIs"<<endl;
			cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1);//npair set to one for now, will be set seperately in execute to the correct value
		}
		else
			cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1 //npair set to one for now, will be set seperately in execute to the correct value
				,WnFolderName,TnFolderName,OutputTnFolderName);

		cosebis->setWns(config->n_max);
		config->cosebis=cosebis;
		return (void *) config;
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

		// clog<<"*************************************"<<endl;
		// clog<<"input_section_name is:"<<config->input_section_name<<endl;
		// clog<<"*************************************"<<endl;

		int num_z_bin_A;
		int num_z_bin_B;

		//get number of tomographic bins
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
		else
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
		}

		vector<number> ell,logell;
		//get ell vector
		status = block->get_val(config->input_section_name, string("ell"), ell);
		if (status) 
		{
			clog<<"Could not load ell in C_ell to COSEBIs"<<endl;
			return status;
		}

		int nell=ell.size();
		//make logell to send to COSEBIs
		for(int i=0; i<nell; i++)
		{
			logell.push_back(log(ell[i]));
		}

		//read input Cls
		vector <vector<number> > InputPower_vec_vec;

		int nPairs=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				// read in C(l)
				vector<number> C_ell;
				string name_in="bin_"+toString(j_bin)+"_"+toString(i_bin);
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
		matrix En_mat=config->cosebis->calEn();
		int n_max=config->n_max;

		if(config->add_cterm)
		{
			number c1=0.;
			number c2=0.;
			block->get_val(cterm_section, "c1", c1);
			block->get_val(cterm_section, "c2", c2);
			matrix En_sin4phi(n_max*nPairs),En_cos4phi(n_max*nPairs);
			int nMaximum_phi=int(config->En_cos4phi.rows/nPairs);
        	for(int i=0;i<n_max*nPairs;i++)
        	{
            	En_cos4phi.load(i,config->En_cos4phi.get((ceil(i/n_max)*nMaximum_phi+i%n_max)));
            	En_sin4phi.load(i,config->En_sin4phi.get((ceil(i/n_max)*nMaximum_phi+i%n_max)));
        	}
			En_mat+=(c1*c1-c2*c2)*En_cos4phi+2.*c1*c2*En_sin4phi;
		}

		if(config->add_2D_cterm)
		{
			number Ac=0.;
			block->get_val(cterm_section, "Ac", Ac);
			matrix En_2D(n_max*nPairs);
			int nMaximum_phi=int(config->En_2D.rows/nPairs);
        	for(int i=0;i<n_max*nPairs;i++)
        	{
            	En_2D.load(i,config->En_2D.get((ceil(i/n_max)*nMaximum_phi+i%n_max)));
        	}
			En_mat+=Ac*En_2D;
		}

  //       //testing m-bias correction
		// int nBins=num_z_bin_A;
		// vector <number> correction(nPairs);
		// if((config->apply_m_bias) && (config->m_bias.size()==nBins))
		// {
		// 	// clog<<"applying m_bias correction, LIKELIHOOD 1"<<endl;
		// 	for(int bin1=0;bin1<nBins;bin1++)
		// 	{
		// 		for(int bin2=bin1;bin2<nBins;bin2++)
		// 		{
		// 			int p1=config->cosebis->calP(nBins,bin1,bin2);
		// 			correction[p1]=(1.+config->m_bias[bin1])*(1.+config->m_bias[bin2]);
		// 			clog<<"p1="<<p1<<" correction="<<correction[p1]<<endl;
		// 			for(int n1=n_max*p1 ;n1<n_max*(p1+1) ;n1++)
		// 			{
		// 				//clog<<"n1="<<n1<<endl;
		// 				En_mat.load(n1,En_mat.get(n1)*correction[p1]);
		// 			}
		// 		}
		// 	}
		// }

		vector<number> En_vec(n_max);
		vector<int> n_vals(n_max);
		int p1=0;
		string name_En;
		for(int i_bin=0; i_bin<num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<num_z_bin_B; j_bin++) 
			{
				int m=0;
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

		if(config->dolikelihood)
		{
			// make a resized data and covariance matrix to match the n_max
	        int nMaximum=int(config->En_data.rows/nPairs);
	        //clog<<"nMaximum="<<nMaximum<<" n_max="<<n_max<<" nPairs"<<nPairs<<endl;
	        if(n_max>nMaximum)
	        {
	            clog<<"set n_max="<<n_max<<" is larger than the input n_max in the data vector:"<<nMaximum<<endl;
	            exit(1);
	        }
	        matrix En_data(n_max*nPairs);
	        //En_data=config->En_data;
	        for(int i=0;i<n_max*nPairs;i++)
	        {
	            En_data.load(i,config->En_data.get((ceil(i/n_max)*nMaximum+i%n_max)));
	        }

	        matrix Cov_mat(n_max*nPairs,n_max*nPairs);
	        //Cov_mat=config->Cov_mat;
	        for(int i=0;i<n_max*nPairs;i++)
				for(int j=0;j<n_max*nPairs;j++)
					Cov_mat.load(i,j,config->Cov_mat.get((ceil(i/n_max)*nMaximum+i%n_max),(ceil(j/n_max)*nMaximum+j%n_max)));

			///Calculate likelihood here
			number ChiS=calChiS(En_mat,En_data,Cov_mat);
			number likelihood_val=-ChiS/2.0;
			string likename=config->output_section_name+"_like";
			status = block->put_val<number>(LIKELIHOODS_SECTION, likename, likelihood_val);
		}
        
		//En_data.printOut((string("En_data_withm")+config->output_section_name+string(".ascii")).c_str(),5);
		//En_mat.printOut((string("En_mat")+config->output_section_name+string(".ascii")).c_str(),5);
		//Cov_mat.printOut((string("Cov_mat")+config->output_section_name+string(".ascii")).c_str(),5);

		return status;
	}

}// end of extern C



    
