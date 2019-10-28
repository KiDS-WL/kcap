///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to xipm and covariance
Under construction
*/
#include "KsiCovariance.h"



extern "C" {
	const string sectionName = "2pcfs";
	const string shear_cl = SHEAR_CL_SECTION;

	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name
	//Some are initialized in the setup, such as xipm. 

	struct pcfs_config {
		string sectionName;
		string input_section_name;
		string input_section_name_plus;
		string input_section_name_minus;
		string output_section_name;
		string Cov_name;

		//read in the min and max theta for the full range of xi_pm
		number theta_min_plus; //minimum theta for xi_plus
		number theta_max_plus; //maximum theta for xi_plus
		number theta_min_minus; //minimum theta for xi_minus
		number theta_max_minus; //maximum theta for xi_minus

		//read the bin centers from file, no bining done for the theory in this case
		matrix theta_mat_plus; //theta values for xi_plus. 
		matrix theta_mat_minus; //theta values for xi_minus. 

		//if weighted binning for the theory is required give edges of each bin.
		vector <number> theta_min_plus_vec; //lower bound for theta_plus per theta bin
		vector <number> theta_max_plus_vec; //upper bound for theta_plus per theta bin
		vector <number> theta_min_minus_vec; //lower bound for theta_plus per theta bin
		vector <number> theta_max_minus_vec; //upper bound for theta_plus per theta bin

		matrix pcfs_data;      //input data 
		matrix Cov_mat;        //input covariance
		//vector<vector<matrix> > weight_plus_mat_vecvec;
		//vector<vector<matrix> > weight_minus_mec_vecvec;
		//string Binning;
		bool constant_cterm_modelling_xim; //if true add cterm modelling
		matrix Sin4phi; // needed for constant cterm modelling for xi_minus
		matrix Cos4phi; // needed for constant cterm modelling for xi_minus
		bool cterm_2D_modelling; //if true add 2D cterm modelling via one scailing parameter Ac
		matrix Xipm_2D_cterm; //Input xipm coming from the 2D cterm 
		bool weight_theta_bins_by_theta; //use weighting in the modelling, 
		//in this case we integrate over the theta values from the given input files
		bool weight_theta_bins_from_input; //use wighting in the modelling,
		//in this case we use a theta weighting and a 100 logbins for each theta bin.
		vector<matrix> theta_Npair_mat_vec; // theta values for the number of galaxy pairs in each redshift bin pair.
		vector<matrix> Npair_mat_vec; // number of galaxy pairs in each redshift bin pair
		vector<matrix> index_min_plus; // we read the input thetas and find the min and max indeces
		vector<matrix> index_max_plus; // given the theta_min_plus/minus_mat and theta_max_plus/minus_mat
		vector<matrix> index_min_minus; // values
		vector<matrix> index_max_minus;

		number sigma_m;
		int nBins;
		bool calCov,calNoiseCov;
		KsiCovariance *xipm;
	} ;
	pcfs_config config;

	

  	void *setup(cosmosis::DataBlock *options) 
  	{
  		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		config.sectionName=OPTION_SECTION;
		string sectionName=config.sectionName;
		
		status=options->get_val<string>(sectionName, string("output_section_name"), string("pcfs_results"), 
			config.output_section_name);
		status=options->get_val<string>(sectionName, string("input_section_name_plus"), "shear_xi_plus", 
			config.input_section_name_plus);
		status=options->get_val<string>(sectionName, string("input_section_name_minus"), "shear_xi_minus", 
			config.input_section_name_minus);

		string theta_plus_file_name;
		string theta_minus_file_name;
		int nTheta_plus,nTheta_minus;
		status=options->get_val<string>(sectionName, string("theta_plus_file_name"), theta_plus_file_name);
		if(status)
		{
			clog<<"Could not load theta_plus_file_name to 2pcfs_interface"<<endl;
			clog<<"looking for theta_min_plus, theta_max_plus and nTheta_plus instead."<<endl;
			status=options->get_val<number>(sectionName, string("theta_min_plus") , config.theta_min_plus);
			if (status) 
			{
				clog<<"Could not load theta_min_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<1.<<endl;
				config.theta_min_plus=1.;
			}
			else
				clog<<"got the value of theta_min_plus="<<config.theta_min_plus<<endl;

			status=options->get_val<number>(sectionName, string("theta_max_plus"), config.theta_max_plus);
			if (status) 
			{
				clog<<"Could not load theta_max_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<100.<<endl;
				config.theta_max_plus=100.;
			}
			else
				clog<<"got the value of theta_max_plus="<<config.theta_max_plus<<endl;

			status=options->get_val<int>(sectionName, string("nTheta_plus") , nTheta_plus);
			if (status) 
			{
				clog<<"Could not load nTheta_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<10<<endl;
				nTheta_plus=10;
			}
			else
				clog<<"got the value of nTheta_plus="<<nTheta_plus<<endl;

			config.theta_mat_plus.resize(nTheta_plus);
			for(int itheta=0; itheta<nTheta_plus; itheta++)
			{
				number theta=exp(log(config.theta_min_plus)+log(config.theta_max_plus/config.theta_min_plus)/(nTheta_plus)*(itheta+0.5));
				config.theta_mat_plus.load(itheta,theta);
			}
		}
		else
		{
			clog<<"got the theta_plus_file_name"<<theta_plus_file_name<<endl;
			config.theta_mat_plus.readFromASCII_marika(theta_plus_file_name.c_str());
		}

		status=options->get_val<string>(sectionName, string("theta_minus_file_name"), theta_minus_file_name);
		if(status)
		{
			clog<<"Could not load theta_minus_file_name to 2pcfs_interface"<<endl;
			clog<<"looking for theta_min_minus, theta_max_minus and nTheta_minus instead."<<endl;
			status=options->get_val<number>(sectionName, string("theta_min_minus") , config.theta_min_minus);
			if (status) 
			{
				clog<<"Could not load theta_min_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<1.<<endl;
				config.theta_min_minus=1.;
			}
			else
				clog<<"got the value of theta_min_minus="<<config.theta_min_minus<<endl;

			status=options->get_val<number>(sectionName, string("theta_max_minus"), config.theta_max_minus);
			if (status) 
			{
				clog<<"Could not load theta_max_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<100.<<endl;
				config.theta_max_minus=100.;
			}
			else
				clog<<"got the value of theta_max_minus="<<config.theta_max_minus<<endl;

			status=options->get_val<int>(sectionName, string("nTheta_minus") , nTheta_minus);
			if (status) 
			{
				clog<<"Could not load nTheta_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<10<<endl;
				nTheta_minus=10;
			}
			else
				clog<<"got the value of nTheta_minus="<<nTheta_minus<<endl;

			config.theta_mat_minus.resize(nTheta_minus);
			for(int itheta=0; itheta<nTheta_minus; itheta++)
			{
				number theta=exp(log(config.theta_min_minus)+log(config.theta_max_minus/config.theta_min_minus)/(nTheta_minus)*(itheta+0.5));
				config.theta_mat_minus.load(itheta,theta);
			}
		}
		else
		{
			clog<<"got the theta_minus_file_name"<<theta_minus_file_name<<endl;
			config.theta_mat_minus.readFromASCII_marika(theta_minus_file_name.c_str());
		}
////
		string theta_min_max_plus_filename;
		vector <number> theta_min_plus_vec;
		vector <number> theta_max_plus_vec;
		status= options->get_val<string>(sectionName, string("theta_min_max_plus_filename"), theta_min_max_plus_filename);
		if(status)
		{
			clog<<"no theta_min_max_plus_filename given"<<endl;
			config.weight_theta_bins_from_input=false;
			config.weight_theta_bins_by_theta=false;
		}
		else
		{
			config.weight_theta_bins_from_input=true;
			config.weight_theta_bins_by_theta=true;
			matrix theta_min_max_plus;
			theta_min_max_plus.readFromASCII_marika(theta_min_max_plus_filename.c_str());
			nTheta_plus=theta_min_max_plus.rows;
			for(int i=0; i<nTheta_plus; i++)
			{
				config.theta_min_plus_vec.push_back(theta_min_max_plus.get(0,i));
				config.theta_max_plus_vec.push_back(theta_min_max_plus.get(1,i));
			}
		}

		string theta_min_max_minus_filename;
		status= options->get_val<string>(sectionName, string("theta_min_max_minus_filename"), theta_min_max_minus_filename);
		if(status)
		{
			clog<<"no theta_min_max_minus_filename given"<<endl;
			config.weight_theta_bins_from_input=false;
			config.weight_theta_bins_by_theta=false;
		}
		else
		{
			config.weight_theta_bins_from_input=true;
			config.weight_theta_bins_by_theta=true;
			matrix theta_min_max_minus;
			theta_min_max_minus.readFromASCII_marika(theta_min_max_minus_filename.c_str());
			nTheta_minus=theta_min_max_minus.rows;
			for(int i=0; i<nTheta_minus; i++)
			{
				config.theta_min_minus_vec.push_back(theta_min_max_minus.get(0,i));
				config.theta_max_minus_vec.push_back(theta_min_max_minus.get(1,i));
			}

		}

		//   initialize 2PCFs
		KsiCovariance *xipm = new KsiCovariance();
		clog<<"going to initialize XIPM"<<endl;

		///change this
		//Need to get an input set of theta_plus and theta_minus to calculate the covariance
		//for. 
		xipm->setParameters(config.theta_mat_plus, config.theta_mat_minus);

		config.Cov_name="XIPMCovarianceFromTheory.ascii";
		status=options->get_val<string>(sectionName, string("cov_name"), config.Cov_name);
		//clog<<"cov name is:"<<config.Cov_En_name<<endl;

		int calculateCov;
		status=options->get_val<int>(sectionName, string("calculateCov"), calculateCov);
		//clog<<"calculate cov:"<<calculateCov<<endl;
		if(status)
		{
			clog<<"calculate cov not set: "<<calculateCov<<endl;
			config.calCov=false;
	  	}
		else
		{
			clog<<"calculateCov set to "<<calculateCov<<endl;
			if(calculateCov)
			{
				config.calCov=true;
				vector<number> sigma_e,ngal_effective;
				number Area;
				status=options->get_val<std::vector<number> >(sectionName, string("sigma_e"),sigma_e);
				if(status)
				{
					clog<<"Didn't find sigma_e values for covariance"<<endl;
			  	}
				else
				{
					clog<<"Found "<<sigma_e.size()<<" sigma_e values"<<endl;
					for(int i=0; i<sigma_e.size(); i++)
					{
						clog<<i<<":"<<sigma_e[i]<<endl;
						sigma_e[i]*=sqrt(2.);
					}
				}
				status=options->get_val<std::vector<number> >(sectionName, string("ngal_effective"),ngal_effective);
				if(status)
				{
					clog<<"Didn't find ngal_effective values for covariance"<<endl;
			  	}
				else
				{
					clog<<"Found "<<ngal_effective.size()<<" ngal_effective values"<<endl;
					for(int i=0; i<ngal_effective.size(); i++)
					{
						clog<<i<<":"<<ngal_effective[i]<<endl;
						ngal_effective[i]*=1./arcmin/arcmin;
					}
				}
				status=options->get_val<number>(sectionName, string("Area"),Area);
				if(status)
				{
					clog<<"Didn't find Area values for covariance"<<endl;
			  	}
				else
				{
					clog<<"Found Area="<<Area<<endl;
					Area*=pow(pi/180.,2);//change to radians
				}
				xipm->setNoise(Area,sigma_e,ngal_effective);
			}
			else
			{
				config.calCov=false;
			}
		}

		config.sigma_m=0.;
		status=options->get_val<number>(sectionName, string("sigma_m"),config.sigma_m);
		if(status)
		{
			clog<<"sigma_m is not set"<<endl;
		}
		else
		{
			clog<<"got the value of sigma_m="<<config.sigma_m<<endl;
		}

		status=options->get_val(sectionName, string("nBins"), config.nBins);
	    if (status) 
	    {
			clog<<"Could not load nBins to XIPM"<<endl;
			// status=get_val("wl_number_density", string("nbin"),config.nBins);
			// clog<<"nBins="<<config.nBins<<endl;
	    }
		else
			clog<<"got the value of nBins="<<config.nBins<<endl;

		string input_nPair_files_suffix;
		status=options->get_val(sectionName, string("input_nPair_files_suffix"),input_nPair_files_suffix);
		if(status)
		{
			config.calNoiseCov=false;
			clog<<"no input_nPair_files_suffix was given"<<endl;
		}
		else
		{
			config.calNoiseCov=true;
			clog<<input_nPair_files_suffix<<" is the input nPair suffix"<<endl;
			if(config.nBins)
			{
				vector<string> FileName_vec;
				for(int bin1=0; bin1<config.nBins; bin1++)
				{
					for(int bin2=bin1; bin2<config.nBins; bin2++)
					{
						string FileName=input_nPair_files_suffix+
							+("_nBins_")+toString(config.nBins)+string("_Bin")
							+toString(bin1+1)+string("_Bin")+toString(bin2+1);
						FileName_vec.push_back(FileName);
					}
				}
				xipm->readNpairs(FileName_vec);
			}
		}
		//exit(1);
		status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_cl"),
			config.input_section_name);
		config.xipm=xipm;
  		return &config;
  		// config is sent to execute 
	}

	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
		// Config is whatever you returned from setup above
		// Block is the collection of parameters and calculations for
		// this set of cosmological parameters
		
		pcfs_config *config= (pcfs_config*) config_in;
		
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//cout<<"config->"
		clog<<"/////////////////\\\\\\\\\\\\\\\\\\\\\\"<<endl;
		clog<<"I'm in execute in 2PCFs_Cl_cosmosis_inteface"<<endl;
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
		// if ((ell[0]>MinPowerCOSEBIs) || (ell[nell-1]<MaxPowerCOSEBIs) || (nell<PowerTableNumberCOSEBIs))
		// {
		// 	clog<<"*************************************"<<endl;
		// 	clog<<"****************WARNING**************"<<endl;
		// 	clog<<"ell range or number of points not sufficient for COSEBIs"<<endl;
		// 	clog<<"your input ell_min="<<ell[0]<<"  ell_max="<<ell[nell-1]<<" n_ell="<<nell<<endl;
		// 	clog<<" For a higher precision at least use this range:"<<endl;
		// 	clog<<"ell_min="<<MinPowerCOSEBIs<<"  ell_max="
		// 		<<MaxPowerCOSEBIs<<" n_ell="<<PowerTableNumberCOSEBIs<<endl;
		// 	clog<<"*************END OF WARNING***********"<<endl;
		// 	clog<<"*************************************"<<endl;
		// }
		
		if (status) 
		{
			clog<<"Could not load ell in C_ell to 2PCFs"<<endl;
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
					clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in C_ell to 2PCFs"<<endl;
					return status;
				}
				InputPower_vec_vec.push_back(C_ell);
				///put Cl in a vector to send to COSEBIs.
				nPairs++;
			}
		}
		clog<<"nPairs="<<nPairs<<endl;
		//status = block->put_val<vector<number> >(config->output_section_name, string("ell.ascii"), ell);
		config->xipm->setZbins(nPairs);
		config->xipm->setPower(logell,InputPower_vec_vec);

		//En_mat=config->xipm->calEn();

		if(config->sigma_m)
		{
			matrix Cov_sigm=config->xipm->calCovForSigma_m(config->sigma_m);
			Cov_sigm.printOut((config->Cov_name+string("_sigma_m_")+toString(config->sigma_m,4)+(".ascii")).c_str(),20);
		}

		if(config->calCov)
		{
			matrix Cov_En=config->xipm->calCov();

			if(config->calNoiseCov)
			{
				if(((config->nBins)*(config->nBins+1))/2.==nPairs)
				{
					Cov_En.printOut((config->Cov_name+string("_TheoryEn.ascii")).c_str(),20);
					matrix CovNoise=config->xipm->calNoiseCov_fromInputNpair();
					CovNoise.printOut((config->Cov_name+string("_NoiseOnly.ascii")).c_str(),20);
					matrix Cov_Bn=config->xipm->calBCov();
					Cov_Bn.printOut((config->Cov_name+string("_TheoryBn.ascii")).c_str(),20);
					config->xipm->setNoiseToZero();
					matrix Cov_cosmicVar=config->xipm->calCov();
					Cov_cosmicVar.printOut((config->Cov_name+string("_TheoryCosmicVar.ascii")).c_str(),20);
					matrix Cov_mixed=Cov_En-Cov_Bn-Cov_cosmicVar;
					Cov_mixed.printOut((config->Cov_name+string("_TheoryMixed.ascii")).c_str(),20);
					Cov_En=Cov_cosmicVar+Cov_mixed+CovNoise;
					Cov_En.printOut((config->Cov_name+string("_NoiseJustForNoise.ascii")).c_str(),20);
					for(int i=0; i<Cov_mixed.rows;i++)
						for(int j=0; j<Cov_mixed.columns;j++)
						{
							if((Cov_Bn.get(i,j)==0) || (CovNoise.get(i,j)/Cov_Bn.get(i,j)<0.))
								Cov_mixed.load(i,j,Cov_mixed.get(i,j));
							else
								Cov_mixed.load(i,j,sqrt(CovNoise.get(i,j)/Cov_Bn.get(i,j))*Cov_mixed.get(i,j));
						}
					Cov_mixed.printOut((config->Cov_name+string("_MixedNoise.ascii")).c_str(),20);
					Cov_En=Cov_cosmicVar+Cov_mixed+CovNoise;
				}
				else
				{
					clog<<"!!!!WARNING!!!! the number of bins for the noise only case does not match the input number of redshift bins"<<endl;
				}
			}
			Cov_En.printOut((config->Cov_name+string(".ascii")).c_str(),20);
			//string nameCovEn="EnCovarianceFromTheory";
			//status=block->put_val<vector<double> >(config->output_section_name, nameCovEn, Cov_En);
		}

		// int n_max=config->n_max;

		// vector<number> En_vec(n_max);
		// vector<int> n_vals(n_max);
		// int p1=0;
		// string name_En;
		// for(int i_bin=0; i_bin<num_z_bin_A; i_bin++) 
		// {
		// 	for (int j_bin=i_bin; j_bin<num_z_bin_B; j_bin++) 
		// 	{
		// 		int m=0;
		// 		for(int n1=n_max*p1,m=0 ;n1<n_max*(p1+1) ;n1++,m++)
		// 			En_vec[m]=En_mat.get(n1);
		// 		name_En=string("bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
		// 		status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
		// 		p1++;
		// 	}
		// }
		// status = block->put_val<int>(config->output_section_name, string("n_mode"), n_max);
		// status = block->put_val<bool>(config->output_section_name, string("b_modes"), config->IsItBmodes);
		// status = block->put_val<double>(config->output_section_name, string("theta_min"), config->theta_min);
		// status = block->put_val<double>(config->output_section_name, string("theta_max"), config->theta_max);

		// for(int n=0;n<n_max;n++)
		// 	n_vals[n]=n+1;

	 //    status = block->put_val<vector<int> >(config->output_section_name, string("cosebis_n"), n_vals);

	    return status;
	}
}// end of extern C


    
