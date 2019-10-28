///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - COSEBIs
*/
#include "COSEBIs.h"


extern "C" 
{
	const string shear_cl = SHEAR_CL_SECTION;
	const number MinPowerCOSEBIs=0.1;
	const number MaxPowerCOSEBIs=1e6;
	const int PowerTableNumberCOSEBIs=200;

	typedef struct COSEBIs_config {
		string sectionName;
		string input_section_name;
		string output_section_name;
		int n_max;
		number theta_min;
		number theta_max;
		int IsItBmodes;
		bool calCov;
		bool calNoiseCov;
		string Cov_En_name;
		int nBins; //this is only needed if the noise only covariance is to be estimated from input nPair
		COSEBIs *cosebis;
		number sigma_m;
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
  	{

  		COSEBIs_config * config = new COSEBIs_config;

  		string sectionName=OPTION_SECTION;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
		
		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);
		if (status) 
		{
			clog<<"Could not load out_section_name to COSEBIs, ";
			clog<<"setting to default: cosebis"<<endl;
			config->output_section_name=string("cosebis");
		}
		else
			clog<<"got the value of output_section_name:"<<config->output_section_name<<endl;

		//get input section name, default= shear_cl
		status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
		if (status) 
		{
			clog<<"Could not load input_section_name to COSEBIs,";
			clog<<" setting to default: shear_cl"<<endl;
			config->input_section_name=string("shear_cl");
		}
		else
			clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;
	

  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config->theta_min);
  		if (status) 
		{
			clog<<"Could not load theta_min to COSEBIs"<<endl;
			clog<<"setting it to the default value of"<<1.<<endl;
		}
		else
			clog<<"got the value of theta_min="<<config->theta_min<<endl;

	    status=options->get_val<number>(sectionName, string("theta_max"),100., config->theta_max);

	    if (status) 
			clog<<"Could not load theta_max to COSEBIs"<<endl;
		else
			clog<<"got the value of theta_max="<<config->theta_max<<endl;


	    status=options->get_val<int>(sectionName, string("n_max"),10, config->n_max);
	    if (status) 
			clog<<"Could not load n_max to COSEBIs"<<endl;
		else
			clog<<"got the value of n_max="<<config->n_max<<endl;


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

		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
		clog<<"going to initialize COSEBIs"<<endl;
		if(defaultFolderName)
		{
			//npair set to one for now, will be set seperately in execute to the correct value
			cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1);
		}
		else
		{
			//npair set to one for now, will be set seperately in execute to the correct value
			cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1 
				,WnFolderName,TnFolderName,OutputTnFolderName);
		}
		
		///read the input COSEBIs file and covarinace file here
		cosebis->setWns(config->n_max);


		config->Cov_En_name="EnCovarianceFromTheory.ascii";
		status=options->get_val<string>(sectionName, string("cov_name"), config->Cov_En_name);
		//clog<<"cov name is:"<<config->Cov_En_name<<endl;

		int calculateCov;
		status=options->get_val<int>(sectionName, string("calculateCov"), calculateCov);
		//clog<<"calculate cov:"<<calculateCov<<endl;
		if(status)
		{
			clog<<"calculate cov not set: "<<calculateCov<<endl;
			config->calCov=false;
	  	}
		else
		{
			clog<<"calculateCov set to "<<calculateCov<<endl;
			if(calculateCov)
			{
				config->calCov=true;
				vector<number> sigma_e,ngal_effective;
				number Area;
				status=options->get_val<std::vector<number> >(sectionName, string("sigma_e"),sigma_e);
				if(status)
				{
					clog<<"Didn't find sigma_e values for covariance"<<endl;
					exit(1);
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
					exit(1);
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
					exit(1);
			  	}
				else
				{
					clog<<"Found Area="<<Area<<endl;
					Area*=pow(pi/180.,2);//change to radians
				}
				cosebis->setNoise(Area,sigma_e,ngal_effective);
			}
			else
			{
				config->calCov=false;
			}
		}

		config->sigma_m=0.;
		status=options->get_val<number>(sectionName, string("sigma_m"),config->sigma_m);
		if(status)
		{
			clog<<"sigma_m is not set"<<endl;
		}
		else
		{
			clog<<"got the value of sigma_m="<<config->sigma_m<<endl;
		}

		status=options->get_val(sectionName, string("nBins"), config->nBins);
	    if (status) 
	    {
			clog<<"Could not load nBins to COSEBIs"<<endl;
			// status=get_val("wl_number_density", string("nbin"),config->nBins);
			// clog<<"nBins="<<config->nBins<<endl;
	    }
		else
			clog<<"got the value of nBins="<<config->nBins<<endl;

		string input_nPair_files_suffix;
		status=options->get_val(sectionName, string("input_nPair_files_suffix"),input_nPair_files_suffix);
		if(status)
		{
			config->calNoiseCov=false;
			clog<<"no input_nPair_files_suffix was given"<<endl;
		}
		else
		{
			config->calNoiseCov=true;
			clog<<input_nPair_files_suffix<<" is the input nPair suffix"<<endl;
			string input_nPair_files_suffix_end="";
			status=options->get_val(sectionName, string("input_nPair_files_suffix_end"),input_nPair_files_suffix_end);
			if(config->nBins)
			{
				vector<string> FileName_vec;
				for(int bin1=0; bin1<config->nBins; bin1++)
				{
					for(int bin2=bin1; bin2<config->nBins; bin2++)
					{
						string FileName=input_nPair_files_suffix+
							+("_nBins_")+toString(config->nBins)+string("_Bin")
							+toString(bin1+1)+string("_Bin")+toString(bin2+1)+input_nPair_files_suffix_end;
						FileName_vec.push_back(FileName);
					}
				}
				cosebis->readNpairs(FileName_vec);
			}
		}

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
				//string name_Cl=string("Cl_bin_")+toString(j_bin)+string("_")+toString(i_bin)+string(".ascii");
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
		clog<<"nPairs="<<nPairs<<endl;
		//status = block->put_val<vector<number> >(config->output_section_name, string("ell.ascii"), ell);
		config->cosebis->setZbins(nPairs);
		config->cosebis->setPower(logell,InputPower_vec_vec);

		En_mat=config->cosebis->calEn();

		if(config->sigma_m)
		{
			matrix Cov_sigm=config->cosebis->calCovForSigma_m(config->sigma_m);
			Cov_sigm.printOut((config->Cov_En_name+string("_sigma_m_")+toString(config->sigma_m,4)+(".ascii")).c_str(),20);
		}

		if(config->calCov)
		{
			matrix Cov_En=config->cosebis->calCov();

			if(config->calNoiseCov)
			{
				if(((config->nBins)*(config->nBins+1))/2.==nPairs)
				{
					Cov_En.printOut((config->Cov_En_name+string("_TheoryEn.ascii")).c_str(),20);
					matrix CovNoise=config->cosebis->calNoiseCov_fromInputNpair();
					CovNoise.printOut((config->Cov_En_name+string("_NoiseOnly.ascii")).c_str(),20);
					matrix Cov_Bn=config->cosebis->calBCov();
					Cov_Bn.printOut((config->Cov_En_name+string("_TheoryBn.ascii")).c_str(),20);
					config->cosebis->setNoiseToZero();
					matrix Cov_cosmicVar=config->cosebis->calCov();
					Cov_cosmicVar.printOut((config->Cov_En_name+string("_TheoryCosmicVar.ascii")).c_str(),20);
					matrix Cov_mixed=Cov_En-Cov_Bn-Cov_cosmicVar;
					Cov_mixed.printOut((config->Cov_En_name+string("_TheoryMixed.ascii")).c_str(),20);
					Cov_En=Cov_cosmicVar+Cov_mixed+CovNoise;
					Cov_En.printOut((config->Cov_En_name+string("_NoiseJustForNoise.ascii")).c_str(),20);
					for(int i=0; i<Cov_mixed.rows;i++)
						for(int j=0; j<Cov_mixed.columns;j++)
						{
							if((Cov_Bn.get(i,j)==0) || (CovNoise.get(i,j)/Cov_Bn.get(i,j)<0.))
								Cov_mixed.load(i,j,Cov_mixed.get(i,j));
							else
								Cov_mixed.load(i,j,sqrt(CovNoise.get(i,j)/Cov_Bn.get(i,j))*Cov_mixed.get(i,j));
						}
					Cov_mixed.printOut((config->Cov_En_name+string("_MixedNoise.ascii")).c_str(),20);
					Cov_En=Cov_cosmicVar+Cov_mixed+CovNoise;
				}
				else
				{
					clog<<"!!!!WARNING!!!! the number of bins for the noise only case does not match the input number of redshift bins"<<endl;
				}
			}
			Cov_En.printOut((config->Cov_En_name+string(".ascii")).c_str(),20);
			//string nameCovEn="EnCovarianceFromTheory";
			//status=block->put_val<vector<double> >(config->output_section_name, nameCovEn, Cov_En);
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
	    return status;
	}
}// end of extern C


    
