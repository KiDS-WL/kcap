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

/*
Things to be tested:
1- c-term modelling
1.1) additive c-terms
1.2) 2D c-terms
*/

///This one calculates things for either xip or xim
extern "C" {
	//section names for inputs from other cosmosis libraries
	const string shear_xi_plus  = "shear_xi_plus";
	const string shear_xi_minus = "shear_xi_minus";
	const string cterm_section  = "shear_c_bias";

	//define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name
	//Some are initialized in the setup.
	typedef struct pcfs_config 
	{
		//input xi_pm section names, either read from the ini file or set to default values
		string input_section_name;//default value: shear_xi_plus if plus or shear_xi_minus if minus

		string output_section_name;//default value: xi_binned, read from ini file in given.

		string type; //plus: xi_+ or minus: xi_-

		//read in the min and max theta for the full range of xi_pm in arcminutes
		number theta_min; //minimum theta for xi, in arcmin
		number theta_max; //maximum theta for xi, in arcmin

		int nTheta; // number of theta bins for xi (assumes log binning)

		//read the bin centers from file, no bining done for the theory in this case
		matrix theta_mat; //theta values for xi_plus in arcmin

		//if weighted binning for the theory is required give edges of each bin.
		vector <number> theta_min_vec; //lower bound for theta per theta bin, in arcmin
		vector <number> theta_max_vec; //upper bound for theta per theta bin, in arcmin

		bool add_c_term; // Set to true if constant c-term modelling is required
		//bool constant_cterm_modelling_xim; //if true add cterm modelling to xi_minus
		matrix sin4phi; // needed for constant cterm modelling for xi_minus
		matrix cos4phi; // needed for constant cterm modelling for xi_minus
		bool inputSinCos_given; //if input files are given in the ini file then will look at the matrices above
		string input_cos4phi_section_name; //if input files not given will look in the block in this section for cos4phi
		string input_sin4phi_section_name; //if input files not given will look in the block in this section for sin4phi

		bool add_2D_cterm; //if true add 2D cterm modelling via one scaling parameter Ac
		matrix Xipm_2D_cterm; //Input xipm coming from the 2D cterm, either xip or xim
		bool inputXi2D_given; // if the input file is given in the ini file then will look at the matrix above
		string input_2D_section_name; //if input files not given will look in the block in this section for xi_2D

		//use weighting in the modelling
		bool weight_theta_bins_by_theta; //use a theta weighting and a 100 logbins for each theta bin.
		bool weight_theta_bins_from_input;//integrate over the theta values from the given input files

		vector<matrix> theta_Npair_mat_vec; //theta values for the number of galaxy pairs in each redshift bin pair in arcmin
		vector<matrix> Npair_mat_vec; //number of galaxy pairs in each redshift bin pair

		//we read the input thetas and find the min and max indices
		//given the theta_min_plus/minus and theta_max_plus/minus values
		vector<matrix> index_min; 
		vector<matrix> index_max; 

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

		//now lets see what is given in the ini file:
		// clog<<endl<<endl;
		// clog<<"*********in xipm_binned interface setup*********"<<endl;
		// clog<<endl;

		
		status=options->get_val<string>(sectionName, string("type"), config->type);
		if(status)
		{
			clog<<"type is not given, going to set to default: plus"<<endl;
			config->type="plus";
		}
		else
		{
			clog<<"got the values of type:"<<config->type<<endl;
			if((config->type=="plus") || (config->type=="minus"))
			{
				clog<<"type value is recognised"<<endl;
			}
			else
			{
				clog<<"please provide a recognised type: plus or minus"<<endl;
				exit(1);
			}
		}

		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);

		if(status)
		{
			config->output_section_name=string("xi_binned_")+config->type;
			clog<<"setting output_section_name to default: "<<config->output_section_name<<endl;
		}
		else
		{
			clog<<"git the value of output_section_name: "<<config->output_section_name<<endl;
		}

		status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
		if(status)
		{
			clog<<"input_section_name is not given, going to set to default value for type:"<<config->type<<endl;
			if(config->type=="plus")
			{
				config->input_section_name=shear_xi_plus;
			}
			else if(config->type=="minus")
			{
				config->input_section_name=shear_xi_minus;
			}
			else
			{
				clog<<"not a recognised type, exiting now ..."<<endl;
				exit(1);
			}
		}
		else
		{
			clog<<"got the values of config->input_section_name:"<<config->input_section_name<<endl;
		}


		//This is a file with just the mid point of the bin given
		string theta_file_name;

		// if status is zero the variable we were looking for is in the ini file 
		// otherwise status return a nonzero value
		bool theta_edges_calculated_from_theta_min_max=false;
		status=options->get_val<string>(sectionName, string("theta_file_name"), theta_file_name);
		if(status)
		{
			clog<<"Could not load theta_file_name to xipm_interface"<<endl;
			clog<<"looking for a list of thetas"<<endl;
			vector<number> theta_vec;
			status=options->get_val<vector<number> >(sectionName, string("theta_list"), theta_vec);
			if(status)
			{
				clog<<"Could not load the theta_list to xipm_interface"<<endl;
				clog<<"looking for theta_min, theta_max and nTheta instead."<<endl;
				status=options->get_val<number>(sectionName, string("theta_min") , config->theta_min);
				if (status) 
				{
					clog<<"Could not load theta_min to xipm_interface"<<endl;
					clog<<"setting it to the default value of:"<<1.<<endl;
					config->theta_min=1.;
				}
				else
					clog<<"got the value of theta_min="<<config->theta_min<<endl;

				status=options->get_val<number>(sectionName, string("theta_max"), config->theta_max);
				if (status) 
				{
					clog<<"Could not load theta_max to xipm_interface"<<endl;
					clog<<"setting it to the default value of:"<<100.<<endl;
					config->theta_max=100.;
				}
				else
					clog<<"got the value of theta_max="<<config->theta_max<<endl;

				status=options->get_val<int>(sectionName, string("nTheta") , config->nTheta);
				if (status) 
				{
					clog<<"Could not load nTheta to xipm_interface"<<endl;
					clog<<"setting it to the default value of:"<<10<<endl;
					config->nTheta=10;
				}
				else
					clog<<"got the value of nTheta="<<config->nTheta<<endl;

				config->theta_mat.resize(config->nTheta);
				for(int itheta=0; itheta<config->nTheta; itheta++)
				{
					number theta=exp(log(config->theta_min)+log(config->theta_max/config->theta_min)/(config->nTheta)*(itheta+0.5));
					config->theta_mat.load(itheta,theta);
				}

				//lets calculate theta_bin edges for weighted binning
				theta_edges_calculated_from_theta_min_max=true;
				clog<<"Setting theta_bin_edges from theta_min theta_max and nTheta. This will be used for wieghted binning"<<endl;
				for(int i=0; i<config->nTheta; i++)
				{
					number theta_low=exp(log(config->theta_min)+log(config->theta_max/config->theta_min)/(config->nTheta)*(i));
					number theta_high=exp(log(config->theta_min)+log(config->theta_max/config->theta_min)/(config->nTheta)*(i+1.0));
					config->theta_min_vec.push_back(theta_low);
					config->theta_max_vec.push_back(theta_high);
					clog<<"theta_low="<<theta_low<<endl;
					clog<<"theta_high="<<theta_high<<endl;
				}
			}
			else
			{
				//now save to theta_mat
				config->nTheta=theta_vec.size();
				clog<<"Found "<<config->nTheta<<" thetas in the list"<<endl;
				config->theta_mat.resize(config->nTheta);
				for(int itheta=0; itheta<config->nTheta; itheta++)
				{
					clog<<"theta_"<<itheta<<"="<<theta_vec[itheta]<<endl;
					config->theta_mat.load(itheta,theta_vec[itheta]);
				}
			}
		}
		else
		{
			clog<<"got the theta_file_name:"<<theta_file_name<<endl;
			config->theta_mat.readFromASCII_marika(theta_file_name.c_str());
		}



		//This should have the same format as the input xipm file
		int add_2D_cterm=0;
		config->inputXi2D_given=false;
		status=options->get_val<int>(sectionName, string("add_2D_cterm"), add_2D_cterm);
		if(add_2D_cterm>0)
		{
			config->add_2D_cterm=true;
			config->input_2D_section_name="xi_2D";
			status=options->get_val<string>(sectionName, string("input_2D_section_name"), config->input_2D_section_name);
			config->input_2D_section_name+string("_")+config->type;
			if(status)
			{
				clog<<"input section name for xi_2D not given setting it to the default value:";
				clog<<config->input_2D_section_name<<endl;
			}
			else
			{
				clog<<"got the input section name for xi_2D: "<<config->input_2D_section_name<<endl;
			}
			
			string InputXipm_2D_cterm;
			status=options->get_val<string>(sectionName, string("InputXipm_2D_cterm"), InputXipm_2D_cterm);
			if(status)
			{
				clog<<"no input Xipm_2D_cterm file was given, will look in the block for input 2D cterm"<<endl;
			}
			else
			{
				clog<<"got the input Xipm_2D_cterm file: "<<InputXipm_2D_cterm<<endl;
				config->Xipm_2D_cterm.readFromASCII_marika((InputXipm_2D_cterm).c_str());
				config->inputXi2D_given=true;
			}

		}
		else
		{
			config->add_2D_cterm=false;
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
		int add_c_term=0;
		config->inputSinCos_given=false;
		status=options->get_val<int>(sectionName, string("add_c_term"), add_c_term);
		if(add_c_term<1)
		{
			clog<<"not using c_term correction"<<endl;
			config->add_c_term=false;
		}
		else
		{
			clog<<"using c_term correction, will look in shear_c_bias section for values"<<endl;
			config->add_c_term=true;

			if(config->type=="minus")
			{
				config->input_cos4phi_section_name="xim_cos4phi";
				status=options->get_val<string>(sectionName, string("input_cos4phi_section_name"), config->input_cos4phi_section_name);
				if(status)
				{
					clog<<"input section name for cos4phi not given setting it to the default value:";
					clog<<config->input_cos4phi_section_name<<endl;
				}
				else
				{
					clog<<"got the input section name for cos4phi: "<<config->input_cos4phi_section_name<<endl;
				}

				config->input_sin4phi_section_name="xim_sin4phi";
				status=options->get_val<string>(sectionName, string("input_sin4phi_section_name"), config->input_sin4phi_section_name);
				if(status)
				{
					clog<<"input section name for sin4phi not given setting it to the default value:";
					clog<<config->input_sin4phi_section_name<<endl;
				}
				else
				{
					clog<<"got the input section name for sin4phi: "<<config->input_sin4phi_section_name<<endl;
				}

				string InputCos4phi;
				status=options->get_val<string>(sectionName, string("InputCos4phi"), InputCos4phi);
				if(status)
				{
					clog<<"no input InputCos4phi file was given, will look in the block"<<endl;
					//config->constant_cterm_modelling_xim=false;
				}
				else
				{
					clog<<"got the input InputCos4phi file: "<<InputCos4phi<<endl;
					config->cos4phi.readFromASCII_marika((InputCos4phi).c_str());
					
					string InputSin4phi;
					status=options->get_val<string>(sectionName, string("InputSin4phi"), InputSin4phi);
					if(status)
					{
						clog<<"no input InputSin4phi file was given, will look in the block"<<endl;
					}
					else
					{
						clog<<"got the input InputSin4phi file: "<<InputSin4phi<<endl;
						config->sin4phi.readFromASCII_marika((InputSin4phi).c_str());
						config->inputSinCos_given=true;
					}
				}
			}
		}

		//Are we going to do a weighted integral over the bins? For this we need to know the edges of the bins
		//If the user doesn't provide this or theta_min theta_max nTheta are not set will not do weighted binning
		string theta_min_max_filename;
		vector <number> theta_min_vec;

		int weighted_binning=1;
		status= options->get_val<int>(sectionName, string("weighted_binning"), weighted_binning);
		if(weighted_binning)
		{
			clog<<"going to do weighted binning"<<endl;
			clog<<"If theta_min_max_filename is not given will use theta_min, theta_max, nTheta to produce the bin edges"<<endl;
			clog<<"If those are not given then will not do weighted binning"<<endl;
			//this should be a file with two columns with column 1 showing theta_min 
			// and column 2 showing theta_max for each theta bin
			// number of theta bins is kept constant for all redshift bin combinations
			status= options->get_val<string>(sectionName, string("theta_min_max_filename"), theta_min_max_filename);
			if(status)
			{
				clog<<"No theta_min_max_filename given"<<endl;
				if(theta_edges_calculated_from_theta_min_max)
				{
					clog<<"theta edges calculated from theta min max"<<endl;
					config->weight_theta_bins_from_input=true;
					config->weight_theta_bins_by_theta=true;

				}
				else
				{
					clog<<"setting weighted integration to false"<<endl;
					config->weight_theta_bins_from_input=false;
					config->weight_theta_bins_by_theta=false;
				}
			}
			else
			{
				config->weight_theta_bins_from_input=true;
				config->weight_theta_bins_by_theta=true;
				matrix theta_min_max;
				theta_min_max.readFromASCII_marika(theta_min_max_filename.c_str());
				if(config->nTheta!=theta_min_max.rows)
				{
					clog<<"the number of rows in theta_min_max_file does not match nTheta, exiting now ..."<<endl;
					exit(1);
				}
				for(int i=0; i<config->nTheta; i++)
				{
					config->theta_min_vec.push_back(theta_min_max.get(0,i));
					config->theta_max_vec.push_back(theta_min_max.get(1,i));
				}

			}

			//This needs to be in this format:
			//files with names starting with {InputNpair}_nBins_${nBins}_Bin${i}_Bin${j}{InputNpair_suffix}
			int n_redshift_pair,Column_theta,Column_Npair,nBins;
			string InputNpair;
			string InputNpair_suffix="";
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
						clog<<"no input Column_theta was given, setting to 1"<<endl;
						Column_theta=1;
					}
					else
					{
						clog<<"got the input Column_theta="<<Column_theta<<endl;
					}
					
					status=options->get_val<int>(sectionName, string("Column_Npair"), Column_Npair);
					if(status)
					{
						clog<<"no input Column_Npair was given, setting to 8"<<endl;
						Column_Npair=8;
					}
					else
					{
						clog<<"got the input Column_Npair="<<Column_Npair<<endl;
					}

					status=options->get_val<int>(sectionName, string("nBins_in"), nBins);
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

					status=options->get_val<string>(sectionName, string("InputNpair_suffix"), InputNpair_suffix);
					if(status)
					{
						clog<<"InputNpair_suffix was not given"<<endl;
					}
					else
					{
						clog<<"got InputNpair_suffix:"<<InputNpair_suffix<<endl;
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
								+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+InputNpair_suffix;
						clog<<"input Napir file name is: "<<Npair_name<<endl;
						FileName_vec.push_back(Npair_name);
					}
				}
				if(CheckFilesExist(FileName_vec))
				{
					//now read in the input npair files and put them in matrices
					for(int r=0; r<n_redshift_pair; r++)
					{
						matrix Npair_mat_in;
						Npair_mat_in.readFromASCII_marika((FileName_vec[r]).c_str());
						matrix theta=Npair_mat_in.getColumn(Column_theta-1);
						vector<number> theta_vec(theta.rows);
						for(int i=0; i<theta.rows; i++)
						{
							theta_vec[i]=theta.get(i);
						}

						matrix Npair_mat=Npair_mat_in.getColumn(Column_Npair-1);

						matrix index_min(config->nTheta);
						matrix index_max(config->nTheta);

						for(int itheta=0; itheta<config->nTheta; itheta++)
						{
							index_min.load(itheta,find_closest_index(theta_vec,config->theta_min_vec[itheta]));
							index_max.load(itheta,find_closest_index(theta_vec,config->theta_max_vec[itheta]));
						}

						config->index_min.push_back(index_min);
						config->index_max.push_back(index_max);
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
		}
		else
		{
			clog<<"Not going to do weighted binning"<<endl;
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

		clog<<endl<<endl;
		clog<<"*********in xipm_binned interface execute*********"<<endl;
		clog<<endl;

		//shear_xi_plus or shear_xi_minus are the the default names of the input section name. 
		//This is another cosmosis library that calculates the theory value of 2pcfs
		int num_z_bin_A;
		int num_z_bin_B;
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			if(status)
			{
				clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
				return status;
			}
			num_z_bin_B = num_z_bin_A;
			//clog<<"num_z_bin_A/B do not exist"<<endl;
		}
		else
		{
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
			if(status)
			{
				clog<<"looked for nbin_b it is not set"<<endl;
				return status;
			}
		}

		//read theta values from shear_xi sections
		vector<number> theta_in;
		status = block->get_val(config->input_section_name, string("theta"), theta_in);
		if (status) 
		{
			clog<<"Could not load theta_in to xipm likelihood"<<endl;
			return status;
		}
		int nTheta_in=theta_in.size();

		// make vectors for logtheta used for interpolation
		vector<number> logtheta(nTheta_in);
		// change units from radians to arcmin
		for (int itheta=0; itheta<nTheta_in; itheta++)
			logtheta[itheta]=log(theta_in[itheta]/arcmin);

		//read in c-term modelling parameters
		number delta_c=0.;
		number c1=0.;
		number c2=0.;
		if(config->add_c_term)
		{
			if(config->type=="plus")
			{
				bool has_val1 = block->has_val(cterm_section, "c1");
				bool has_val2 = block->has_val(cterm_section, "c2");
				if(has_val1 & has_val2)
				{
					block->get_val(cterm_section, "c1", c1);
					block->get_val(cterm_section, "c2", c2);
					delta_c=c1*c1+c2*c2;
				}
				else
				{
					status=block->get_val(cterm_section, "delta_c", delta_c);
					if(status)
					{
						clog<<"no input c_term is given"<<endl;
						return status;
					}
				}
			}
			else if (config->type=="minus")
			{
				status = block->get_val(cterm_section, "c1", c1);
				if(status)
				{
					clog<<"input c1 is not given in the value.ini, going to exit now..."<<endl;
					return status;
				}
				status = block->get_val(cterm_section, "c2", c2);
				if(status)
				{
					clog<<"input c2 is not given in the value.ini, going to exit now..."<<endl;
					return status;
				}
			}
		}

		number Ac=0.;
		if(config->add_2D_cterm)
		{	
			status=block->get_val(cterm_section, "Ac", Ac);
			if(status)
			{
				clog<<"input Ac is not given, so going to exit ..."<<endl;
				return status;
			}
		}

		///Here is where the binning is done and c-terms are added
		int pair=0;	
		int start_ind=0;
		int start_ind_2D=0;
		
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=1; j_bin<=num_z_bin_B; j_bin++) 
			{
				// read in input xi
				vector<number> xi;
				string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
				string theta_name_in=string("theta_")+name_in;
				
				bool has_val = block->has_val(config->input_section_name, name_in);
				if (has_val) 
				{
					status = block->get_val<vector<number> >(config->input_section_name, name_in, xi);
					function_cosebis pcfs_table;
					pcfs_table.loadWithValues(logtheta,xi,true);
					pcfs_table.extrapolationOff();
					matrix xi_binned(config->nTheta);

					if(config->weight_theta_bins_from_input || config->weight_theta_bins_by_theta)
					{
						//clog<<"going to do weighted binning"<<endl;
						for(int itheta=0; itheta<config->nTheta; itheta++)
						{
							number weightedTheory=0.;
							if (config->weight_theta_bins_by_theta)
							{
								matrix bin_mat(100),weight_mat(100);
								for(int it=0; it<100;it++)
								{
									number tmin=config->theta_min_vec[itheta];
									number tmax=config->theta_max_vec[itheta];
									number theta=exp(log(tmin)+log(tmax/tmin)/(100)*(it+0.5));
									number deltaTheta=exp(log(tmin)+log(tmax/tmin)/(100)*(it+1.0))-exp(log(tmin)+log(tmax/tmin)/(100)*(it+0.));
									bin_mat.load(it,pcfs_table.value(theta));
									weight_mat.load(it,theta*deltaTheta);
								}
								matrix mult=(bin_mat.t()*weight_mat);
								//clog<<"mult size is"<<mult.size();
								weightedTheory=mult.get(0,0)/weight_mat.sum();
							}
							else
							{
								int itmin=int(config->index_min[pair].get(itheta));
								int itmax=int(config->index_max[pair].get(itheta));
								matrix bin_mat(itmax-itmin+1),weight_mat(itmax-itmin+1);
								for(int it=itmin; it <= itmax; it++)
								{
									number theta=config->theta_Npair_mat_vec[pair].get(it);
									bin_mat.load(it-itmin,pcfs_table.value(theta));
									weight_mat.load(it-itmin,config->Npair_mat_vec[pair].get(it));
								}
								matrix mult=(bin_mat.t()*weight_mat);
								//clog<<"mult size is"<<mult.size();
								weightedTheory=mult.get(0,0)/weight_mat.sum();
							}
							xi_binned.load(0,itheta,weightedTheory);
						}
					}// end of weighted binning
					//no weighted binning just use the single theta values for the theory vector
					else
					{
						
						for(int itheta=0; itheta<config->nTheta; itheta++)
							xi_binned.load(0,itheta, pcfs_table.value(config->theta_mat.get(itheta)));
					}
					// we have calculated xi_binned

					// now add c-term contributions
					if(config->add_c_term)
					{
						clog<<"adding cosntant c-term"<<endl;
						if(config->type=="plus")
						{
							matrix ones_mat(config->nTheta);
							matrix xi_plus_c1c2=delta_c*ones_mat;
							xi_binned+=xi_plus_c1c2;
						}
						if(config->type=="minus")
						{
							matrix sin4phi(config->nTheta),cos4phi(config->nTheta);
							if(config->inputSinCos_given)
							{
								for(int ind=0; ind<config->nTheta; ind++)
							 	{
							 		cos4phi.load(ind,config->cos4phi.get(ind+start_ind));
							 		sin4phi.load(ind,config->sin4phi.get(ind+start_ind));
							 	}
							 	start_ind+=config->nTheta;
							}
							else
							{
								bool has_val1 = block->has_val(config->input_cos4phi_section_name, name_in);
								bool has_val2 = block->has_val(config->input_sin4phi_section_name, name_in);
								vector<number> cos4phi_vec, sin4phi_vec;
								if (has_val1 & has_val2) 
								{
									status = block->get_val<vector<number> >(config->input_cos4phi_section_name, name_in, cos4phi_vec);
									status = block->get_val<vector<number> >(config->input_sin4phi_section_name, name_in, sin4phi_vec);
									//now check if the size is correct
									if((cos4phi_vec.size()==config->nTheta) &(sin4phi_vec.size()==config->nTheta))
									{
										for(int i=0; i<config->nTheta; i++)
										{
											cos4phi.load(i,cos4phi_vec[i]);
											sin4phi.load(i,sin4phi_vec[i]);
										}
									}
								}
								else
								{
									clog<<"****WARNING no value for cos4phi or sin4phi given in the block, exiting now ..."<<endl;
									exit(1);
								}
							}
							matrix xi_minus_c1c2=(c1*c1-c2*c2)*cos4phi+2.*c1*c2*sin4phi;	
							xi_binned+=xi_minus_c1c2;
						}
					} //end of add_c_term

					if(config->add_2D_cterm)
					{
						clog<<"adding 2D cterm"<<endl;
						matrix xi_2D(config->nTheta);
						if(config->inputXi2D_given)
						{
							for(int ind=0; ind<config->nTheta; ind++)
						 	{
						 		xi_2D.load(ind,config->Xipm_2D_cterm.get(ind+start_ind_2D));
						 	}
				        	start_ind_2D+=config->nTheta;
						}
						else
						{
							bool has_val = block->has_val(config->input_2D_section_name, name_in);
							vector<number> xi_2D_vec;
							if (has_val) 
							{
								status = block->get_val<vector<number> >(config->input_2D_section_name, name_in, xi_2D_vec);
								//now check if the size is correct
								if(xi_2D_vec.size()==config->nTheta)
								{
									for(int i=0; i<config->nTheta; i++)
									{
										xi_2D.load(i,xi_2D_vec[i]);
									}
								}
							}
							else
							{
								clog<<"*****WARNING no value for xi_2D is given in the block exiting now..."<<endl;
								exit(1);
							}
				        }
				        //Sould this be Ac or Ac^2
						xi_binned+=Ac*xi_2D;
					} //end of add 2D cterm

					//could add the cuts here if need be
					
					vector<number> xi_vec(xi_binned.rows),theta_vec(config->theta_mat.rows);
					for(int m=0 ;m<xi_binned.rows ;m++)
					{
						xi_vec[m]=xi_binned.get(m);
						theta_vec[m]=config->theta_mat.get(m);
					}


					status = block->put_val<vector<number> >(config->output_section_name, name_in, xi_vec);
					status = block->put_val<vector<number> >(config->output_section_name, theta_name_in, theta_vec);
					
					pair++;
				} // end of hasval for input xi
			} //end of second for loop
		} //end of first for loop
		return status;
	} //end of execute
}// end of extern C

    
