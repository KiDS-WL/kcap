///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from xi_pm to binned xi_pm with added effects.
Also can read the data files, cut scales and put them in the same format. 
*/
#include "matrix.h"
#include "function_cosebis.h"
#include "calChiS.h"


extern "C" {
	//section names for inputs from other cosmosis libraries
	const string shear_xi_plus = "shear_xi_plus";
	const string shear_xi_minus = "shear_xi_minus";

		//define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name
	//Some are initialized in the setup.
	typedef struct pcfs_config {
		
		string sectionName;//name of this module given in cosmosis ini file

		//input xi_pm section names, either read from the ini file or set to default values
		string input_section_name_plus;//default value: shear_xi_plus
		string input_section_name_minus;//default value: shear_xi_minus

		string output_section_name;//default value: pcfs_results, read from ini file in given.

		//read in the min and max theta for the full range of xi_pm
		number theta_min_plus; //minimum theta for xi_plus, in arcmin
		number theta_max_plus; //maximum theta for xi_plus, in arcmin
		number theta_min_minus; //minimum theta for xi_minus, in arcmin
		number theta_max_minus; //maximum theta for xi_minus, in arcmin

		int nTheta_minus;
		int nTheta_plus;

		//read the bin centers from file, no bining done for the theory in this case
		matrix theta_mat_plus; //theta values for xi_plus. in arcmin
		matrix theta_mat_minus; //theta values for xi_minus. in arcmin

		//if weighted binning for the theory is required give edges of each bin.
		vector <number> theta_min_plus_vec; //lower bound for theta_plus per theta bin, in arcmin
		vector <number> theta_max_plus_vec; //upper bound for theta_plus per theta bin, in arcmin
		vector <number> theta_min_minus_vec; //lower bound for theta_minus per theta bin, in arcmin
		vector <number> theta_max_minus_vec; //upper bound for theta_minus per theta bin, in arcmin

		matrix pcfs_data;      //input data file
		matrix Cov_mat;        //input covariance

		
		bool constant_cterm_modelling;
		bool constant_cterm_modelling_xim; //if true add cterm modelling to xi_minus
		matrix Sin4phi; // needed for constant cterm modelling for xi_minus
		matrix Cos4phi; // needed for constant cterm modelling for xi_minus

		bool cterm_2D_modelling; //if true add 2D cterm modelling via one scaling parameter Ac
		matrix Xipm_2D_cterm; //Input xipm coming from the 2D cterm 

		//use weighting in the modelling
		bool weight_theta_bins_by_theta; //use a theta weighting and a 100 logbins for each theta bin.
		bool weight_theta_bins_from_input;//integrate over the theta values from the given input files

		vector<matrix> theta_Npair_mat_vec; //theta values for the number of galaxy pairs in each redshift bin pair. in arcmin
		vector<matrix> Npair_mat_vec; //number of galaxy pairs in each redshift bin pair

		//we read the input thetas and find the min and max indices
		//given the theta_min_plus/minus and theta_max_plus/minus values
		vector<matrix> index_min_plus; 
		vector<matrix> index_max_plus; 
		vector<matrix> index_min_minus; 
		vector<matrix> index_max_minus;

		// calculate likelihood if data vector and covariance are given
		bool dolikelihood;

	} pcfs_config; //make an the structure to be used


	// this is needed to find the name of the module in the cosmosis ini file.
	int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
	{
		auto status = options->get_val(OPTION_SECTION, name, parameter);
		if (status!=DBS_SUCCESS) 
		{
			parameter = "";
			cerr<< "Could not find or understand parameter in 2PCFs section: " 
				<< name << std::endl; 
			return 1;
		}
		return 0;
	}

	// this function simply finds the redshift bin pair given nBins 
	// and the first and second redshift bins: fbin and sbin
	// Could change this given a z1 and z2 that are at least 1 and nBins 
	// to return index=nBins*(z1-1)-(z1-1)*(z1-1-1)/2+(z2-1)-(z1-1)
	// but ain't nobody got time for that
	int calP(int nBins,int fbin,int sbin)
	{
		if(fbin>sbin)
		{
			int swap=fbin;
			fbin=sbin;
			sbin=swap;
		}
		int p=fbin*nBins;

		for(int i=0; i<fbin; i++)
		  	p-=i;
		return p+sbin-fbin;
	}

	// here is the setup function. Cosmosis runs this only once and sends the address to a structure with
	// all that is needed to execute. 
	void *setup(cosmosis::DataBlock *options) 
	{
		//options reads the ini file
		//define config here and then read from options the relevant input quantities
  		pcfs_config * config = new pcfs_config;

		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//set the section name from the ini file
		string sectionName=OPTION_SECTION;
		config->sectionName=sectionName;
		//now lets see what is given in the ini file

		//get_val takes values from the ini file in the given sectionName.
		//get_val<string>(sectionName, "name of the variable to read", "some default value for it", VariableToBeSet);
		status=options->get_val<string>(sectionName, string("output_section_name"), string("pcfs"), 
			config->output_section_name);
		status=options->get_val<string>(sectionName, string("input_section_name_plus"), shear_xi_plus, 
			config->input_section_name_plus);
		status=options->get_val<string>(sectionName, string("input_section_name_minus"), shear_xi_minus, 
			config->input_section_name_minus);

		string theta_plus_file_name;
		string theta_minus_file_name;
		//int nTheta_plus,nTheta_minus;
		// if status is zero the variable we were looking for is in the ini file 
		// otherwise status return a nonzero value
		status=options->get_val<string>(sectionName, string("theta_plus_file_name"), theta_plus_file_name);
		if(status)
		{
			clog<<"Could not load theta_plus_file_name to 2pcfs_interface"<<endl;
			clog<<"looking for theta_min_plus, theta_max_plus and nTheta_plus instead."<<endl;
			status=options->get_val<number>(sectionName, string("theta_min_plus") , config->theta_min_plus);
			if (status) 
			{
				clog<<"Could not load theta_min_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<1.<<endl;
				config->theta_min_plus=1.;
			}
			else
				clog<<"got the value of theta_min_plus="<<config->theta_min_plus<<endl;

			status=options->get_val<number>(sectionName, string("theta_max_plus"), config->theta_max_plus);
			if (status) 
			{
				clog<<"Could not load theta_max_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<100.<<endl;
				config->theta_max_plus=100.;
			}
			else
				clog<<"got the value of theta_max_plus="<<config->theta_max_plus<<endl;

			status=options->get_val<int>(sectionName, string("nTheta_plus") , config->nTheta_plus);
			if (status) 
			{
				clog<<"Could not load nTheta_plus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<10<<endl;
				config->nTheta_plus=10;
			}
			else
				clog<<"got the value of nTheta_plus="<<config->nTheta_plus<<endl;

			config->theta_mat_plus.resize(config->nTheta_plus);
			for(int itheta=0; itheta<config->nTheta_plus; itheta++)
			{
				number theta=exp(log(config->theta_min_plus)+log(config->theta_max_plus/config->theta_min_plus)/(config->nTheta_plus)*(itheta+0.5));
				config->theta_mat_plus.load(itheta,theta);
			}
		}
		else
		{
			clog<<"got the theta_plus_file_name"<<theta_plus_file_name<<endl;
			config->theta_mat_plus.readFromASCII_marika(theta_plus_file_name.c_str());
		}

		status=options->get_val<string>(sectionName, string("theta_minus_file_name"), theta_minus_file_name);
		if(status)
		{
			clog<<"Could not load theta_minus_file_name to 2pcfs_interface"<<endl;
			clog<<"looking for theta_min_minus, theta_max_minus and nTheta_minus instead."<<endl;
			status=options->get_val<number>(sectionName, string("theta_min_minus") , config->theta_min_minus);
			if (status) 
			{
				clog<<"Could not load theta_min_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<1.<<endl;
				config->theta_min_minus=1.;
			}
			else
				clog<<"got the value of theta_min_minus="<<config->theta_min_minus<<endl;

			status=options->get_val<number>(sectionName, string("theta_max_minus"), config->theta_max_minus);
			if (status) 
			{
				clog<<"Could not load theta_max_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<100.<<endl;
				config->theta_max_minus=100.;
			}
			else
				clog<<"got the value of theta_max_minus="<<config->theta_max_minus<<endl;

			status=options->get_val<int>(sectionName, string("nTheta_minus") , config->nTheta_minus);
			if (status) 
			{
				clog<<"Could not load nTheta_minus to 2pcfs_interface"<<endl;
				clog<<"setting it to the default value of:"<<10<<endl;
				config->nTheta_minus=10;
			}
			else
				clog<<"got the value of nTheta_minus="<<config->nTheta_minus<<endl;

			config->theta_mat_minus.resize(config->nTheta_minus);
			for(int itheta=0; itheta<config->nTheta_minus; itheta++)
			{
				number theta=exp(log(config->theta_min_minus)+log(config->theta_max_minus/config->theta_min_minus)/(config->nTheta_minus)*(itheta+0.5));
				config->theta_mat_minus.load(itheta,theta);
			}
		}
		else
		{
			clog<<"got the theta_minus_file_name"<<theta_minus_file_name<<endl;
			config->theta_mat_minus.readFromASCII_marika(theta_minus_file_name.c_str());
		}

		//Are we going to do a weighted integral over the bins? For this we need to know the edges of the bins
		//If the user doesn't provide this we don't integrate
		string theta_min_max_plus_filename;
		string theta_min_max_minus_filename;
		vector <number> theta_min_plus_vec;
		vector <number> theta_max_plus_vec;
		//this should be a file with two columns with column 1 showing theta_min 
		//and column two showing theta_max for each theta bin
		// number of theta bins is kept constant for all redshift bin combinations
		// if you want different number of bins for each redshift bin pair then you need to
		// modify the code
		status= options->get_val<string>(sectionName, string("theta_min_max_plus_filename"), theta_min_max_plus_filename);
		if(status)
		{
			clog<<"no theta_min_max_plus_filename given, setting weighted integration to false"<<endl;
			config->weight_theta_bins_from_input=false;
			config->weight_theta_bins_by_theta=false;
		}
		else
		{
			config->weight_theta_bins_from_input=true;
			config->weight_theta_bins_by_theta=true;
			matrix theta_min_max_plus;
			theta_min_max_plus.readFromASCII_marika(theta_min_max_plus_filename.c_str());
			config->nTheta_plus=theta_min_max_plus.rows;
			for(int i=0; i<config->nTheta_plus; i++)
			{
				config->theta_min_plus_vec.push_back(theta_min_max_plus.get(0,i));
				config->theta_max_plus_vec.push_back(theta_min_max_plus.get(1,i));
			}

			//now check it for xi_m
			status= options->get_val<string>(sectionName, string("theta_min_max_minus_filename"), theta_min_max_minus_filename);
			if(status)
			{
				clog<<"no theta_min_max_minus_filename given, setting weighted integration to false"<<endl;
				config->weight_theta_bins_from_input=false;
				config->weight_theta_bins_by_theta=false;
			}
			else
			{
				config->weight_theta_bins_from_input=true;
				config->weight_theta_bins_by_theta=true;
				matrix theta_min_max_minus;
				theta_min_max_minus.readFromASCII_marika(theta_min_max_minus_filename.c_str());
				config->nTheta_minus=theta_min_max_minus.rows;
				for(int i=0; i<config->nTheta_minus; i++)
				{
					config->theta_min_minus_vec.push_back(theta_min_max_minus.get(0,i));
					config->theta_max_minus_vec.push_back(theta_min_max_minus.get(1,i));
				}
			}
		}


      	string Input2PCFsFileName;
		status=options->get_val<string>(sectionName, string("input_2pcfs_filename"), Input2PCFsFileName);
		if(status)
		{
			clog<<"no input 2PCFs file was given, not going to calculate likelihoods"<<endl;
			config->dolikelihood=false;
		}
		else
		{
			clog<<"got the input 2PCFs files: "<<Input2PCFsFileName<<endl;
			config->pcfs_data.readFromASCII_marika((Input2PCFsFileName).c_str());
			config->dolikelihood=true;

			///read Covariance from file which should correspond to xi_p ,xi_m for each redshift pair 
			///The data vector format should be:
			// \xi_+^{11}(\theta_1)
			// \xi_+^{11}(\theta_2)
			// ...
			// \xi_+^{11}(\theta_max)
			// \xi_+^{nn}(\theta_1)
			// ...
			// \xi_+^{nn}(\theta_max)
			// \xi_-^{11}(\theta_1)
			// \xi_-^{11}(\theta_2)
			// ...
			// \xi_-^{11}(\theta_max)
			// ...
			// \xi_-^{nn}(\theta_1)
			// ...
			// \xi_-^{nn}(\theta_max)
		    string InputCovarianceFileName;
			status=options->get_val<string>(sectionName, string("input_covariance_filename"), InputCovarianceFileName);
			if(status)
			{
				clog<<"no input Covariance file was given, not going to calculate likelihood"<<endl;
				config->dolikelihood=false;
			}
			else
			{
				clog<<"got the input covariance file: "<<InputCovarianceFileName<<endl;
				config->Cov_mat.readFromASCII_marika((InputCovarianceFileName).c_str());
				config->dolikelihood=true;
			}
		}

		//This should have the same format as the input xipm file
		string InputXipm_2D_cterm;
		status=options->get_val<string>(sectionName, string("InputXipm_2D_cterm"), InputXipm_2D_cterm);
		if(status)
		{
			clog<<"no input Xipm_2D_cterm file was given"<<endl;
			config->cterm_2D_modelling=false;
		}
		else
		{
			clog<<"got the input Xipm_2D_cterm file: "<<InputXipm_2D_cterm<<endl;
			config->Xipm_2D_cterm.readFromASCII_marika((InputXipm_2D_cterm).c_str());
			config->cterm_2D_modelling=true;
		}

		//this should be in this format:
		// cos4phi^{11}(\theta_1)
		// cos4phi^{11}(\theta_2)
		// ...
		// cos4phi^{11}(\theta_max)
		// ...
		// cos4phi^{nn}(\theta_1)
		// ...
		// cos4phi^{nn}(\theta_max)
		int add_c_term;
		status=options->get_val<int>(sectionName, string("add_c_term"), add_c_term);
		if(status || (add_c_term<1))
		{
			clog<<"not using c_term correction"<<endl;
			config->constant_cterm_modelling=false;
		}
		else
		{
			clog<<"using c_term correction, will look in shear_c_bias section for values"<<endl;
			config->constant_cterm_modelling=true;

			string InputCos4phi;
			status=options->get_val<string>(sectionName, string("InputCos4phi"), InputCos4phi);
			if(status)
			{
				clog<<"no input InputCos4phi file was given"<<endl;
				config->constant_cterm_modelling_xim=false;
			}
			else
			{
				clog<<"got the input InputCos4phi file: "<<InputCos4phi<<endl;
				config->Cos4phi.readFromASCII_marika((InputCos4phi).c_str());
				config->constant_cterm_modelling_xim=true;
				string InputSin4phi;
				status=options->get_val<string>(sectionName, string("InputSin4phi"), InputSin4phi);
				if(status)
				{
					clog<<"no input InputSin4phi file was given"<<endl;
					config->constant_cterm_modelling_xim=false;
				}
				else
				{
					clog<<"got the input InputSin4phi file: "<<InputSin4phi<<endl;
					config->Sin4phi.readFromASCII_marika((InputSin4phi).c_str());
					config->constant_cterm_modelling_xim=true;
				}
			}
		}

		

		//This needs to be in this format:
		//files with names starting with {InputNpair}_nBins_${nBins}_Bin${i}_Bin${j}
		int n_redshift_pair,Column_theta,Column_Npair,nBins;
		string InputNpair;
		if(config->weight_theta_bins_from_input)
		{
			status=options->get_val<string>(sectionName, string("InputNpair"), InputNpair);
			if(status)
			{
				clog<<"no input InputNpair file was given, going to set weighting to theta"<<endl;
				config->weight_theta_bins_from_input=false;
			}
			else
			{
				clog<<"got the input InputNpair file: "<<InputNpair<<endl;

				status=options->get_val<int>(sectionName, string("Column_theta"), Column_theta);
				if(status)
				{
					clog<<"no input Column_theta was given, setting to 0"<<endl;
					Column_theta=0;
				}
				else
				{
					clog<<"got the input Column_theta="<<Column_theta<<endl;
				}
				
				status=options->get_val<int>(sectionName, string("Column_Npair"), Column_Npair);
				if(status)
				{
					clog<<"no input Column_Npair was given, setting to 7"<<endl;
					Column_Npair=7;
				}
				else
				{
					clog<<"got the input Column_Npair="<<Column_Npair<<endl;
				}

				status=options->get_val<int>(sectionName, string("nBins"), nBins);
				if(status)
				{
					clog<<"no input nBins was given, can't read input Npair files"<<endl;	
					config->weight_theta_bins_from_input=false;
				}
				else
				{
					clog<<"got the number of redshift bins="<<nBins<<endl;
					n_redshift_pair=int((nBins+1.)*(nBins)/2.);
					config->weight_theta_bins_from_input=true;
				}
			}
		}

		if(config->weight_theta_bins_from_input)
		{
			// set the name of the files to be read
			vector<string> FileName_vec;
			for (int i_bin=1; i_bin<=nBins; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=nBins; j_bin++) 
				{
					string Npair_name="";
					Npair_name=InputNpair
							+string("_nBins_")+toString(nBins)
							+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin);
					clog<<"input Napir file name is: "<<Npair_name<<endl;
					FileName_vec.push_back(Npair_name);
				}
			}
			if(CheckFilesExist(FileName_vec))
			{
				//now read in the input napir files and put them in matrices
				for(int r=0; r<n_redshift_pair; r++)
				{
					matrix Npair_mat_in;
					Npair_mat_in.readFromASCII_marika((FileName_vec[r]).c_str());
					matrix theta=Npair_mat_in.getColumn(Column_theta);
					vector<number> theta_vec(theta.rows);
					for(int i=0; i<theta.rows; i++)
					{
						theta_vec[i]=theta.get(i);
					}
					matrix Npair_mat=Npair_mat_in.getColumn(Column_Npair);

					matrix index_min_plus(config->nTheta_plus);
					matrix index_max_plus(config->nTheta_plus);

					matrix index_min_minus(config->nTheta_minus);
					matrix index_max_minus(config->nTheta_minus);

					for(int itheta=0; itheta<config->nTheta_plus; itheta++)
					{
						index_min_plus.load(itheta,find_closest_index(theta_vec,config->theta_min_plus_vec[itheta]));
						index_max_plus.load(itheta,find_closest_index(theta_vec,config->theta_max_plus_vec[itheta]));
					}
					for(int itheta=0; itheta<config->nTheta_minus; itheta++)
					{
						index_min_minus.load(itheta,find_closest_index(theta_vec,config->theta_min_minus_vec[itheta]));
						index_max_minus.load(itheta,find_closest_index(theta_vec,config->theta_max_minus_vec[itheta]));
					}
					config->index_min_plus.push_back(index_min_plus);
					config->index_max_plus.push_back(index_max_plus);
					config->index_min_minus.push_back(index_min_minus);
					config->index_max_minus.push_back(index_max_minus);
					config->theta_Npair_mat_vec.push_back(theta);
					config->Npair_mat_vec.push_back(Npair_mat);
				}
				config->weight_theta_bins_by_theta=false;
			}
			else
			{
				clog<<"one or more Npair files don't exist, goin to use Npair proportional to theta instead"<<endl;
				config->weight_theta_bins_from_input=false;
				config->weight_theta_bins_by_theta=true;
			}
		}
		
		return (void *) config;
		// config is sent to execute 
	}

	//This is where the calculations are made. Everything that is setup and sent through config is used here
	// This runs at every step of the sampler
	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
		// Config is whatever you returned from setup above
		// Block is the collection of parameters and calculations for
		// this set of cosmological parameters
		pcfs_config *config= (pcfs_config*) config_in;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//For now only works for nbin_a=nbin_b
		int num_z_bin_A;
		int num_z_bin_B;
		//shear_xi_plus is the the default name of the input section name. 
		//This is another cosmosis library that calculates the theory value of 2pcfs
		status = block->get_val(config->input_section_name_plus, string("nbin_a"), num_z_bin_A);
		if(status)
		{
			status = block->get_val(config->input_section_name_plus, string("nbin_b"), num_z_bin_B);
		}
		else
		{
			status = block->get_val(config->input_section_name_plus, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
		}


		vector<number> theta_plus_in,theta_minus_in;
		status = block->get_val(config->input_section_name_plus, string("theta"), theta_plus_in);
		if (status) 
		{
			clog<<"Could not load theta_plus to 2PCFS liklihood"<<endl;
			return status;
		}
		int nTheta_plus_in=theta_plus_in.size();

		status = block->get_val(config->input_section_name_minus, string("theta"), theta_minus_in);	
		if (status) 
		{
			clog<<"Could not load theta_minus to 2PCFS liklihood"<<endl;
			return status;
		}
		int nTheta_minus_in=theta_minus_in.size();

		//These are vectors of vectors, effectively they work as matrices.
		vector <vector<number> > Input_xip_vec_vec;
		vector <vector<number> > Input_xim_vec_vec;
		int nPairs=0;

		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				// read in 2pcfs
				vector<number> xip;
				vector<number> xim;
				//these used to be different now they are the same. 
				string name_in_plus="bin_"+toString(j_bin)+"_"+toString(i_bin);
				string name_in_minus="bin_"+toString(j_bin)+"_"+toString(i_bin);
	    		status = block->get_val(config->input_section_name_plus, name_in_plus, xip);
	    		if (status) 
				{
					clog<<"Could not load xip bin "<<j_bin<<"_"<< i_bin<<" in 2PCFS likelihood"<<endl;
					return status;
				}

	    		status = block->get_val(config->input_section_name_minus, name_in_minus, xim);
				if (status) 
				{
					clog<<"Could not load xim bin "<<j_bin<<"_"<< i_bin<<" in 2PCFS likelihood"<<endl;
					return status;
				}
				//push the new vectors to the end of the double vectors
				Input_xip_vec_vec.push_back(xip);
				Input_xim_vec_vec.push_back(xim);
				nPairs++;
			}
		}
		
		///Make a vector of xip and xim based on the theta cuts 
		// and integrate over theta for the binning as well.
		int nTheta_plus,nTheta_minus;
		// If weighted binning is true 
		// then the number of theta bins is set by theta_min_plus_vec and theta_min_minus_vec
		if(config->weight_theta_bins_from_input || config->weight_theta_bins_by_theta)
		{
			nTheta_plus  =config->theta_min_plus_vec.size();
			nTheta_minus =config->theta_min_minus_vec.size();
		}
		// If weighted binning is false 
		// then number of theta bins is set by theta_mat_plus and theta_mat_minus
		else
		{
			nTheta_plus  =config->theta_mat_plus.size();
			nTheta_minus =config->theta_mat_minus.size();
		}
		
		//make a matrix object for the theory vector
		matrix pcfs_mat((nTheta_minus+nTheta_plus)*nPairs);
		// make vectors for logtheta used for interpolation
		vector<number> logtheta_minus(nTheta_minus_in),logtheta_plus(nTheta_plus_in);
		// change units from radians to arcmin
		for (int itheta=0; itheta<nTheta_minus_in; itheta++)
		{
			logtheta_minus[itheta]=log(theta_minus_in[itheta]/arcmin);
		}

		for (int itheta=0; itheta<nTheta_plus_in; itheta++)
			logtheta_plus[itheta]=log(theta_plus_in[itheta]/arcmin);

		// vector<number> theta_vec_plus(nTheta_plus);
		// for (int itheta=0; itheta<nTheta_plus; itheta++)
		// 	theta_vec_plus[itheta]=config->theta_mat_plus.get(itheta);

		// vector<number> theta_vec_minus(nTheta_minus);
		// for (int itheta=0; itheta<nTheta_minus; itheta++)
		// 	theta_vec_minus[itheta]=config->theta_mat_minus.get(itheta);

		
		int counterP=0;
		int redshift=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				//makes a function object. This can tabulate values and interpolate them among other things
				function_cosebis pcfs_table;
				vector<number> xi_p_vec(nTheta_plus);

				pcfs_table.loadWithValues(logtheta_plus,Input_xip_vec_vec[redshift],true);
				pcfs_table.extrapolationOff();

				//do weigthed binning if true
				if(config->weight_theta_bins_from_input || config->weight_theta_bins_by_theta)
				{
					for(int itheta=0; itheta<nTheta_plus; itheta++)
					{
						number weightedTheory=0.;
						if (config->weight_theta_bins_by_theta)
						{
							matrix bin_mat(100),weight_mat(100);
							for(int it=0; it<100;it++)
							{
								number theta=exp(log(config->theta_min_plus_vec[itheta])+log(config->theta_max_plus_vec[itheta]/config->theta_min_plus_vec[itheta])/(100)*(it+0.5));
								bin_mat.load(it,pcfs_table.value(theta));
								weight_mat.load(it,theta);
							}
							matrix mult=(bin_mat.t()*weight_mat);
							//clog<<"mult size is"<<mult.size();
							weightedTheory=mult.get(0,0)/weight_mat.sum();
						}
						else
						{
							// vector<matrix> theta_Npair_mat_vec; // theta values for the number of galaxy pairs in each redshift bin pair.
							// vector<matrix> Npair_mat_vec; // number of galaxy pairs in each redshift bin pair
							// vector<matrix> index_min_plus; // we read the input thetas and find the min and max indeces
							// vector<matrix> index_max_plus; // given the theta_min_plus/minus_mat and theta_max_plus/minus_mat
							// vector<matrix> index_min_minus; // values
							// vector<matrix> index_max_minus;
							int itmin=int(config->index_min_plus[redshift].get(itheta));
							int itmax=int(config->index_max_plus[redshift].get(itheta));
							matrix bin_mat(itmax-itmin+1),weight_mat(itmax-itmin+1);
							for(int it=itmin; it <= itmax; it++)
							{
								number theta=config->theta_Npair_mat_vec[redshift].get(it);
								bin_mat.load(it-itmin,pcfs_table.value(theta));
								weight_mat.load(it-itmin,config->Npair_mat_vec[redshift].get(it));
							}
							matrix mult=(bin_mat.t()*weight_mat);
							//clog<<"mult size is"<<mult.size();
							weightedTheory=mult.get(0,0)/weight_mat.sum();
						}
						pcfs_mat.load(0,counterP,weightedTheory);
						xi_p_vec[itheta]=pcfs_mat.get(0,counterP);
						counterP++;
					}
				}
				//no weighted binning just use the single theta values for the theory vector
				else
				{
					for(int itheta=0; itheta<nTheta_plus; itheta++)
					{
						pcfs_mat.load(0,counterP, pcfs_table.value(config->theta_mat_plus.get(itheta)));
						xi_p_vec[itheta]=pcfs_mat.get(0,counterP);
						counterP++;
					}
				}
				//to save or not? Have to check this doesn't happen all the time.
				string name_xi_p=string("xi_p_bin_")+toString(j_bin)+string("_")+toString(i_bin);
				status = block->put_val<vector<double> >(config->output_section_name, name_xi_p, xi_p_vec);
				redshift++;
			}
		}
		
		//do the same for xim
		redshift=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				function_cosebis pcfs_table;
				vector<number> xi_m_vec(nTheta_minus);
				pcfs_table.loadWithValues(logtheta_minus,Input_xim_vec_vec[redshift],true);
				pcfs_table.extrapolationOff();

				if(config->weight_theta_bins_from_input || config->weight_theta_bins_by_theta)
				{
					for(int itheta=0; itheta<nTheta_minus; itheta++)
					{
						number weightedTheory=0.;
						if (config->weight_theta_bins_by_theta)
						{
							matrix bin_mat(100),weight_mat(100);
							for(int it=0; it<100;it++)
							{
								number theta=exp(log(config->theta_min_minus_vec[itheta])+log(config->theta_max_minus_vec[itheta]/config->theta_min_minus_vec[itheta])/(100)*(it+0.5));
								bin_mat.load(it,pcfs_table.value(theta));
								weight_mat.load(it,theta);
							}
							matrix mult=(bin_mat.t()*weight_mat);
							//clog<<"mult size is"<<mult.size();
							weightedTheory=mult.get(0,0)/weight_mat.sum();
						}
						else
						{
							// vector<matrix> theta_Npair_mat_vec; // theta values for the number of galaxy pairs in each redshift bin pair.
							// vector<matrix> Npair_mat_vec; // number of galaxy pairs in each redshift bin pair
							// vector<matrix> index_min_plus; // we read the input thetas and find the min and max indeces
							// vector<matrix> index_max_plus; // given the theta_min_plus/minus_mat and theta_max_plus/minus_mat
							// vector<matrix> index_min_minus; // values
							// vector<matrix> index_max_minus;
							int itmin=int(config->index_min_minus[redshift].get(itheta));
							int itmax=int(config->index_max_minus[redshift].get(itheta));
							matrix bin_mat(itmax-itmin+1),weight_mat(itmax-itmin+1);
							for(int it=itmin; it <= itmax; it++)
							{
								number theta=config->theta_Npair_mat_vec[redshift].get(it);
								bin_mat.load(it-itmin,pcfs_table.value(theta));
								weight_mat.load(it-itmin,config->Npair_mat_vec[redshift].get(it));
							}
							matrix mult=(bin_mat.t()*weight_mat);
							//clog<<"mult size is"<<mult.size();
							weightedTheory=mult.get(0,0)/weight_mat.sum();
						}
						pcfs_mat.load(0,counterP,weightedTheory);
						xi_m_vec[itheta]=pcfs_mat.get(0,counterP);
						counterP++;
					}
				}
				else
				{
					for(int itheta=0; itheta<nTheta_minus; itheta++)
					{
						pcfs_mat.load(0,counterP, pcfs_table.value(config->theta_mat_minus.get(itheta)));
						xi_m_vec[itheta]=pcfs_mat.get(0,counterP);
						counterP++;
					}
				}
				//to save or not
				string name_xi_m=string("xi_m_bin_")+toString(j_bin)+string("_")+toString(i_bin);
				status = block->put_val<vector<double> >(config->output_section_name, name_xi_m, xi_m_vec);
				redshift++;
			}
		}
		
		//2D cterm
		if(config->cterm_2D_modelling)
		{
			number Ac=0.;
			status = block->get_val<number>("shear_c_bias", "Ac", Ac);
			pcfs_mat+=Ac*config->Xipm_2D_cterm;
		}

		//add cterm modelling for both xim and xip
		if(config->constant_cterm_modelling)
		{
			number c1=0.0;
			number c2=0.0;
			status = block->get_val<number>("shear_c_bias", "c1", c1);
			if(status)
			{
				//clog<<"no c1 in shear_c_bias"<<endl;
				status = block->get_val<number>("shear_c_bias", "c", c1);
				if(status)
				{
					clog<<"no c1 or c foound in shear_c_bias"<<endl;
					return status;
				}
				else
				{
					//clog<<"setting c1=c2"<<endl;
					c2=c1;
					config->constant_cterm_modelling_xim=false;
				}
			}
			else
			{
				clog<<"looking for c2"<<endl;
				status = block->get_val<number>("shear_c_bias", "c2", c2);
				if(status)
				{
					clog<<"didn't find c2"<<endl;
					return status;
				}
			}
			matrix ones_mat(nTheta_plus*nPairs);
			matrix xi_plus_c1c2=(c1*c1+c2*c2)*ones_mat;
			matrix xi_c1c2((nTheta_minus+nTheta_plus)*nPairs);
			for (int i=0; i<nTheta_plus*nPairs; i++)
				xi_c1c2.load(i,xi_plus_c1c2.get(i));

			//add cterm modelling for both xim and xip
			if(config->constant_cterm_modelling_xim)
			{
				matrix xi_minus_c1c2=(c1*c1-c2*c2)*config->Cos4phi+2.*c1*c2*config->Sin4phi;
				for (int i=0; i<nTheta_minus*nPairs; i++)
					xi_c1c2.load(i+nTheta_plus*nPairs,xi_minus_c1c2.get(i));
			}
			else
			{
				for (int i=0; i<nTheta_minus*nPairs; i++)
					xi_c1c2.load(i+nTheta_plus*nPairs,0);
			}
			pcfs_mat+=xi_c1c2;
		}




		

		//now save xipm to block, have to be written
		// vector<number> xip(nTheta_plus);
		// vector<number> xim(nTheta_minus);
		// int p1=0;
		// string name;
		// for(int i_bin=0; i_bin<num_z_bin_A; i_bin++) 
		// {
		// 	for (int j_bin=i_bin; j_bin<num_z_bin_B; j_bin++) 
		// 	{
		// 		int m=0;
		// 		for(int n1=nTheta_plus*p1,m=0 ;n1<nTheta_plus*(p1+1) ;n1++,m++)
		// 			xip_vec[m]=xi_p_vec[n1];

		// 		m=0;
		// 		for(int n1=nTheta_minus*p1,m=0 ;n1<(nTheta_minus*(p1+1)) ;n1++,m++)
		// 			xim_vec[m]=xi_m_vec.[n1];

		// 		name=string("bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
		// 		status = block->put_val<vector<double> >(config->output_section_name+string("_plus"), name, xip);
		// 		status = block->put_val<vector<double> >(config->output_section_name+string("_minus"), name, xim);
		// 		p1++;
		// 	}
		// }
		
		status = block->put_val<vector<double> >(config->output_section_name+string("_plus"), string("theta_min_plus_vec"), config->theta_min_plus_vec);
		status = block->put_val<vector<double> >(config->output_section_name+string("_plus"), string("theta_max_plus_vec"), config->theta_max_plus_vec);
		status = block->put_val<vector<double> >(config->output_section_name+string("_minus"), string("theta_min_minus_vec"), config->theta_min_minus_vec);
		status = block->put_val<vector<double> >(config->output_section_name+string("_minus"), string("theta_max_minus_vec"), config->theta_max_minus_vec);

		status = block->put_val<double>(config->output_section_name+string("_plus"), string("theta_min_plus"), config->theta_min_plus);
		status = block->put_val<double>(config->output_section_name+string("_plus"), string("theta_max_plus"), config->theta_max_plus);
		status = block->put_val<int>(config->output_section_name+string("_plus"), string("nTheta_plus"), config->nTheta_plus);

		status = block->put_val<double>(config->output_section_name+string("_minus"), string("theta_min_minus"), config->theta_min_minus);
		status = block->put_val<double>(config->output_section_name+string("_minus"), string("theta_max_minus"), config->theta_max_minus);
		status = block->put_val<int>(config->output_section_name+string("_minus"), string("nTheta_minus"), config->nTheta_minus);

		//These few lines actually calculate the likelihood
		if(config->dolikelihood)
		{
			// matrix Delta=config->pcfs_data-pcfs_mat;
			// Delta.printOut(string("Delta").c_str(),5);
			number ChiS=calChiS(pcfs_mat,config->pcfs_data,config->Cov_mat);
			//number ChiS=chiS.get(0);
			//clog<<"Delta.rows="<<Delta.rows<<endl;
			//clog<<"Delta.cols="<<Delta.columns<<endl;
			number likelihood_val=-ChiS/2.0;
			string likename=config->output_section_name+"_like";
			status = block->put_val<number>(LIKELIHOODS_SECTION, likename, likelihood_val);
		}
		return status;
	}
}// end of extern C

    
