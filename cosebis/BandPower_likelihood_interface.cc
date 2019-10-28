///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - BandPower
*/
#include "BandPower.h"
#include "calChiS.h"


extern "C" 
{

	typedef struct BandPower_config {
		string sectionName;
		string input_section_name;
		string output_section_name;
		number thetamin;
		number thetamax;
		vector<vector<int> >index_vec_vec; //this has the indices for the bins to be calculated for each z1 z2 combination given in z1_z2_mat
		matrix z1_z2_mat;// this has all the z1 z2 combinations that should be calculated
		//matrix cutsmat; //this matrix has the cuts for each redshift bin saved in it. It is read from an input file
		bool docuts; //if cutsmat is given then do cuts
		vector<number> l_min_vec;
		vector<number> l_max_vec;
		int IsItBmodes;
		int type; //0: clustering, 2: GGL, 3: cosmic shear E-modes, 4: cosmic shear B-modes
		//bool calCov;
		//bool calNoiseCov;
		//string Cov_En_name;
		int nBands;
		BandPower *BP0;
		BandPower *BP2;
		BandPower *BP4;
		number sigma_m;
		matrix BP_data;
		matrix Cov_mat;
		bool dolikelihood;
	}
	BandPower_config;

	//const string shear_cl = SHEAR_CL_SECTION;
	const int NLBINS=5000;
	const number LLOW=0.1;
	const number LHIGH= 1e4;

	//define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as BandPower. 
	//BandPower is a class that produces P_E and P_B for cosmic shear, ggl and clustering using different 
	//weight functions. Covariance is not caluclated here yet. 

	//find out what is the name of this section in the ini file. 
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

	//setup only runs once when cosmosis is called
  	void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
  	{
		BandPower_config * config = new BandPower_config;

		string sectionName=OPTION_SECTION;
		config->sectionName=sectionName;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
		
		// This is where the outputs will be saved
		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);
		if (status) 
		{
			clog<<"Could not load out_section_name to BandPower, ";
			clog<<"setting to default: bandpower"<<endl;
			config->output_section_name=string("bandpower");
		}
		else
			clog<<"got the value of output_section_name:"<<config->output_section_name<<endl;

		// This is where the input data file is. If data or covariance are not given likelihood is not caluclated
		string input_bandpower_filename;
		status=options->get_val<string>(sectionName, string("input_bandpower_filename"), input_bandpower_filename);
		if (status) 
		{
			clog<<"Could not load input_bandpower_filename to BandPower, ";
			clog<<"will not produce likelihood"<<endl;
			config->dolikelihood=false;
		}
		else
		{
			clog<<"got the value of input_bandpower_filename:"<<input_bandpower_filename<<endl;
			///read BP_data from file
			config->BP_data.readFromASCII_marika((input_bandpower_filename).c_str());
			config->dolikelihood=true;
			//get input covariance file name
			//if covariance is not given likelihood is not calculated
			string InputCovarianceFileName;
			status=options->get_val<string>(sectionName, string("input_covariance_filename"), InputCovarianceFileName);
			if (status) 
			{
				clog<<"Could not find the name of the input Covariance for bandpower file"<<endl;
				config->dolikelihood=false;
			}
			else
			{
				clog<<"got the name of the input Covariance for BandPower file="<<InputCovarianceFileName<<endl;
				///read Covariance from file
	        	config->Cov_mat.readFromASCII_marika((InputCovarianceFileName).c_str());
	        }
		}
	
		//minimum theta value that is used to estimate band powers from xi_pm, in arcmins
		//make sure this matches what is used for the data
  		status=options->get_val<number>(sectionName, string("theta_min") , config->thetamin);
  		if (status) 
		{
			clog<<"Could not load theta_min to BandPower"<<endl;
			clog<<"setting it to the default value of 1 arcmins"<<endl;
			config->thetamin=1.;
		}
		else
			clog<<"got the value of theta_min="<<config->thetamin<<endl;

		//maximum theta value that is used to estimate band powers from xi_pm, in arcmins
		//make sure this matches what is used for the data.
	    status=options->get_val<number>(sectionName, string("theta_max"), config->thetamax);
	    if (status) 
	    {
			clog<<"Could not load theta_max to BandPower, setting it to the default value of 100 arcmin"<<endl;
			config->thetamax=100.;
	    }
		else
			clog<<"got the value of theta_max="<<config->thetamax<<endl;

		//change them to randian
		config->thetamin*=arcmin;
		config->thetamax*=arcmin;

		string cuts_file;
		status=options->get_val<string>(sectionName, string("cuts_per_redshift_bin"), cuts_file);
		if(status)
		{
			config->docuts=false;
			clog<<"Could not load cuts_per_redshift_bin to BandPower"; 
		}
		else
		{
			clog<<"reading cuts_file"<<endl;
			matrix cutsmat;
			cutsmat.readFromASCII_marika(cuts_file.c_str());
			config->docuts=true;
			vector<int> rows_to_remove;
			vector<int> columns_to_remove;
			for(int col=2; col<cutsmat.columns; col++)
			{
				columns_to_remove.push_back(col);
				clog<<"columns to remove:"<<col<<endl;
			}
			clog<<cutsmat.rows<<"  "<<cutsmat.columns<<endl;
			//now make a matrix with all z1 z2 combinations in the cutcov.

			int z1=0;
			int z2=0;
			int pair=0;
			vector<int> index_vec;
			for (int i_cut=0; i_cut<cutsmat.rows; i_cut++)
			{
				clog<<"i_cut="<<i_cut<<endl;
				int z1_temp=cutsmat.get(0,i_cut);
				int z2_temp=cutsmat.get(1,i_cut);
				clog<<"z1_temp="<<z1_temp<<" z2_temp="<<z2_temp<<endl;
				int index=cutsmat.get(2,i_cut)-1;
				if((z1_temp==z1) && (z2_temp==z2))
				{
					clog<<"zi and zi_temp match"<<endl;
					rows_to_remove.push_back(i_cut);
					clog<<"rows to remove:"<<i_cut<<endl;
					index_vec.push_back(index);
				}
				else
				{
					if(i_cut>0)
					{
						pair++;
						config->index_vec_vec.push_back(index_vec);
						// matrix index_mat(index_vec.size());
						// for(int a=0; a<index_vec.size(); a++)
						// 	index_mat.load(a,index_vec[a]);
						//index_mat.printOut((string("index_mat")+toString(pair)).c_str(),0);
						index_vec.clear();
					}
					index_vec.push_back(index);
					z1=z1_temp;
					z2=z2_temp;
				}
			}
			pair++;
			config->index_vec_vec.push_back(index_vec);
			// matrix index_mat(index_vec.size());
			// for(int a=0; a<index_vec.size(); a++)
			// {
			// 	index_mat.load(a,index_vec[a]);
			// }
			//index_mat.printOut((string("index_mat")+toString(pair)).c_str(),0);
			index_vec.clear();
			config->z1_z2_mat=subMatrix_removeRows_and_Columns(cutsmat,rows_to_remove,columns_to_remove);
			//z1_z2_mat.printOut("z1_z2_mat.ascii",0);

		}
		//exit(1);

		//first lets look for a file with l_min l_max for each band power bin
		string l_min_max_file;
		//vectors that contain the min and max value of ell for each band power bin
		//vector<number> l_min_vec,l_max_vec;
		status=options->get_val<string>(sectionName, string("l_min_max_file"), l_min_max_file);
		//int nBands;
		if (status)
		{
			number l_min,l_max;
			clog<<"Could not load l_min_max_file to BandPower"; 
			clog<<"going to look for nBands and l_min and l_max instead"<<endl;
			//number of bands
		    status=options->get_val<int>(sectionName, string("nBands"), config->nBands);
		    if (status) 
				clog<<"Could not load n_max to BandPower, setting it to the default value of 10"<<endl;
			else
				clog<<"got the value of nBands="<<config->nBands<<endl;

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
			
			for(int b=0; b<config->nBands; b++)
			{
				number lmin=exp(log(l_min)+log(l_max/l_min)/(config->nBands)*b);
				number lmax=exp(log(l_min)+log(l_max/l_min)/(config->nBands)*(b+1));
				config->l_min_vec.push_back(lmin);
				config->l_max_vec.push_back(lmax);
				clog<<"b="<<b<<" lmin="<<lmin<<" lmax="<<lmax<<endl;
			}
		}
		else
		{
			matrix l_min_max_mat;
			l_min_max_mat.readFromASCII_marika(l_min_max_file.c_str());
			config->nBands=l_min_max_mat.rows;
			for(int b=0; b<config->nBands; b++)
			{
				config->l_min_vec.push_back(l_min_max_mat.get(0,b));
				config->l_max_vec.push_back(l_min_max_mat.get(1,b));
			}
		}
		
		//If set to 1 then calculates B-mode band powers, needs b-mode Cls as input.
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


		//this is where the bandpower weights will be saved
		string FolderName;
		status=options->get_val<string>(sectionName, string("Output_FolderName"), FolderName);
		if(status)
		{
			FolderName="./BandPower_weights";
			clog<<"Could not find the folder name for weight functions in W_Output_FolderName, setting to default:";
			clog<<FolderName<<endl;
		}
		else
			clog<<"folder name for bandpower weights is:"<<FolderName<<endl;

		//is set for apodisation.
		number  Delta_x;
		status=options->get_val<number>(sectionName, string("Delta_x"), Delta_x);
		if(status)
		{
			clog<<"Could not find Delta_x, setting to default:"<<endl;
			Delta_x=0.5;
		}
		else
			clog<<"Delta_theta is:"<<Delta_x<<endl;

		//This is the bandpower filter type that we want to mimic.
		//The only option that the code understands for now is tophat
		string Response_function_type;
		status=options->get_val<string>(sectionName, string("Response_function_type"), Response_function_type);
		if(status)
		{
			clog<<"Could not find the response function type setting to default: tophat"<<endl;
			Response_function_type="tophat";
		}
		else
			clog<<"response function type is:"<<Response_function_type<<endl;

		// use analytic solution for ?
		bool Analytic=true; //set this to be read from the ini file
		int Analytic_int;
		status=options->get_val<int>(sectionName, string("Analytic"), Analytic_int);
		if(status)
		{
			clog<<"Could not find Analytic, setting to default: true"<<endl;
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

		//Should I apodise?
		bool noApodise=true; //set this to be read from the ini file
		int Apodise;
		string ap_str="";
		status=options->get_val<int>(sectionName, string("Apodise"), Apodise);
		if(status)
		{
			clog<<"Could not find Apodise, setting to default: false"<<endl;
			ap_str="_noAp";
		}
		else
		{
			if(Apodise==0)
			{
				noApodise=true;
				clog<<"setting apodise to false"<<endl;
				ap_str="_noAp";
			}
			else
			{
				noApodise=false;
				clog<<"setting apodise to true"<<endl;
				ap_str="_ap";
			}
		}
		
		//set file names for g and W filters
		string gFileName="g";
		string WFileName="W"+ap_str;
		string type_str;
		bool real=false;
		//   initialize BandPower based on the type the user asks for
		// options are: clustering, ggl, cosmic_shear_e and cosmic_shear_b
		status=options->get_val<string>(sectionName, string("type"), type_str);
		if(status)
		{
			clog<<"Could not find the band power type to calculate, setting to default"<<endl;
			config->type=0;
		}
		else
		{
			clog<<"band power type to calculate is:"<<type_str<<endl;
			if(type_str=="clustering")
			{
				config->type=0;
				BandPower *BP0 = new BandPower();
				BP0->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						0,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);
				config->BP0=BP0;

				//check if input section name is given in the ini file otherwise use the default
				status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
				if (status) 
				{
					clog<<"Could not load input_section_name to band power,";
					clog<<" setting to default: "<< GALAXY_CL_SECTION <<endl;
					config->input_section_name=GALAXY_CL_SECTION;
				}
				else
					clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;
			}
			else if(type_str=="ggl")
			{
				config->type=2;
				int bessel_order=2;
				BandPower *BP2 = new BandPower();
				BP2->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						bessel_order,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);
				config->BP2=BP2;
				//check if input section name is given in the ini file otherwise use the default
				status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
				if (status) 
				{
					clog<<"Could not load input_section_name to band power,";
					clog<<" setting to default: "<< "galaxy_shear_cl" <<endl;
					config->input_section_name="galaxy_shear_cl";
				}
				else
					clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;
			}
			else if(type_str=="cosmic_shear_e")
			{
				config->type=3;
				BandPower *BP0 = new BandPower();
				BP0->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						0,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);

				config->BP0=BP0;
				BandPower *BP4 = new BandPower();
				BP4->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						4,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);
				config->BP4=BP4;
				//check if input section name is given in the ini file otherwise use the default
				status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
				if (status) 
				{
					clog<<"Could not load input_section_name to band power,";
					clog<<" setting to default: "<< SHEAR_CL_SECTION <<endl;
					config->input_section_name=SHEAR_CL_SECTION;
				}
				else
					clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;

			}
			else if(type_str=="cosmic_shear_b")
			{
				config->type=4;
				BandPower *BP0 = new BandPower();
				BP0->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						0,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);
				config->BP0=BP0;

				BandPower *BP4 = new BandPower();
				BP4->initialize(config->thetamin,config->thetamax,Response_function_type,
						config->l_min_vec,config->l_max_vec,
						4,noApodise,Delta_x, Analytic,
						LLOW,
						LHIGH,
						NLBINS,
						FolderName,gFileName,WFileName,real);
				config->BP4=BP4;
				//****Note: B-mode CL is not included in cosmosis. For now this is set to E-mode CL***
				//check if input section name is given in the ini file otherwise use the default
				status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
				if (status) 
				{
					clog<<"Could not load input_section_name to band power,";
					clog<<" setting to default: "<< SHEAR_CL_SECTION <<endl;
					config->input_section_name=SHEAR_CL_SECTION;
				}
				else
					clog<<"got the value of input_section_name:"<<config->input_section_name<<endl;
			}
			else
			{
				clog<<"not a recognised type, please select from: clustering, ggl, cosmic_shear_e or cosmic_shear_b"<<endl;
				clog<<"exiting now ..."<<endl;
				exit(1);
			}
		}
		
		clog<<"BP initialised ....."<<endl;

		return (void *) config;
  		//return &config;
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
		int num_z_bin_A;
		int num_z_bin_B;
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
			clog<<"num_z_bin_A/B do not exist"<<endl;
			
		}
		else
		{
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
		}

		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"num_z_bin_A="<<num_z_bin_A<<" num_z_bin_B="<<num_z_bin_B<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		// clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		

		vector<number> ell,logell;
		status = block->get_val(config->input_section_name, string("ell"), ell);
		if (status) 
		{
			clog<<"Could not load ell in "<<config->input_section_name<<" to BandPower"<<endl;
			return status;
		}

		int nell=ell.size();
		for(int i=0; i<nell; i++)
			logell.push_back(log(ell[i]));


		if(config->docuts)
		{
			for(int pair=0; pair<config->z1_z2_mat.rows; pair++)
			{
				int z1=config->z1_z2_mat.get(0,pair);
				int z2=config->z1_z2_mat.get(1,pair);
				vector<number> C_ell;
				string name_in=string("bin_")+toString(z1)+string("_")+toString(z2);
				string index_name_in=string("index_")+name_in;
				bool has_val = block->has_val(config->input_section_name, name_in);
				if (has_val) 
				{
					status = block->get_val<vector<number> >(config->input_section_name, name_in, C_ell);
					//find BP for this bin and save it.
					matrix BP_mat;
					//clog<<"setting BP"<<endl;
					if(config->type==0)
					{
						config->BP0->setInput_single(logell,C_ell);
						BP_mat=config->BP0->calBP(config->index_vec_vec[pair]);
						clog<<"BP_mat.rows="<<BP_mat.rows<<endl;
					}
					else if(config->type==2)
					{
						config->BP2->setInput_single(logell,C_ell);
						BP_mat=config->BP2->calBP(config->index_vec_vec[pair]);
					}
					else if(config->type==3)
					{
						config->BP0->setInput_single(logell,C_ell);
						matrix BP_mat0;
						BP_mat0=config->BP0->calBP(config->index_vec_vec[pair]);

						config->BP4->setInput_single(logell,C_ell);
						matrix BP_mat4;
						BP_mat4=config->BP4->calBP(config->index_vec_vec[pair]);

						BP_mat=BP_mat0/2.+BP_mat4/2.;
					}
					else if(config->type==4)
					{
						config->BP0->setInput_single(logell,C_ell);
						matrix BP_mat0;
						BP_mat0=config->BP0->calBP(config->index_vec_vec[pair]);

						config->BP4->setInput_single(logell,C_ell);
						matrix BP_mat4;
						BP_mat4=config->BP4->calBP(config->index_vec_vec[pair]);

						BP_mat=BP_mat0/2.-BP_mat4/2.;
					}
					vector<number> BP_vec(BP_mat.rows);
					for(int m=0 ;m<BP_mat.rows ;m++)
						BP_vec[m]=BP_mat.get(m);
					status = block->put_val<vector<number> >(config->output_section_name, name_in, BP_vec);
					status = block->put_val<vector<int> >(config->output_section_name, index_name_in, config->index_vec_vec[pair]);
				}
				else
				{
					clog<<"Could not load bin_"<<z1<<"_"<<z2<<" in C_ell to bandpower, exiting now"<<endl;
					exit(1);
				}
			}
		}// if config->docuts end
		else
		{
			vector<int> n_bands(config->nBands);
			for(int n=0;n<config->nBands;n++)
			 	n_bands[n]=n+1;

			for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
			{
				for (int j_bin=1; j_bin<=num_z_bin_B; j_bin++) 
				{
					// read in C(l)
					vector<number> C_ell;
					string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
					string index_name_in=string("index_")+name_in;
					bool has_val = block->has_val(config->input_section_name, name_in);
					//clog<<"name_in="<<name_in<<endl;
					if (has_val) 
					{
						status = block->get_val<vector<number> >(config->input_section_name, name_in, C_ell);
						//find BP for this bin and save it.
						matrix BP_mat;
						//clog<<"setting BP"<<endl;
						if(config->type==0)
						{
							config->BP0->setInput_single(logell,C_ell);
							BP_mat=config->BP0->calBP();
						}
						else if(config->type==2)
						{
							config->BP2->setInput_single(logell,C_ell);
							BP_mat=config->BP2->calBP();
						}
						else if(config->type==3)
						{
							config->BP0->setInput_single(logell,C_ell);
							matrix BP_mat0;
							BP_mat0=config->BP0->calBP();

							config->BP4->setInput_single(logell,C_ell);
							matrix BP_mat4;
							BP_mat4=config->BP4->calBP();

							BP_mat=BP_mat0/2.+BP_mat4/2.;
						}
						else if(config->type==4)
						{
							config->BP0->setInput_single(logell,C_ell);
							matrix BP_mat0;
							BP_mat0=config->BP0->calBP();

							config->BP4->setInput_single(logell,C_ell);
							matrix BP_mat4;
							BP_mat4=config->BP4->calBP();

							BP_mat=BP_mat0/2.-BP_mat4/2.;
						}
						number nBands=BP_mat.rows;
						vector<number> BP_vec(nBands);
						for(int m=0 ;m<nBands ;m++)
							BP_vec[m]=BP_mat.get(m);
						status = block->put_val<vector<number> >(config->output_section_name, name_in, BP_vec);
						status = block->put_val<vector<int> >(config->output_section_name, index_name_in, n_bands);
					}
					else
					{
						clog<<"has_val is false. Could not load bin_"<<j_bin<<"_"<< i_bin<<" in C_ell to bandpower"<<endl;
						//return status;
					}
				}
			}
		}


		status = block->put_val<vector<number> >(config->output_section_name, string("l_min_vec"), config->l_min_vec);
		status = block->put_val<vector<number> >(config->output_section_name, string("l_max_vec"), config->l_max_vec);
		status = block->put_val<bool>(config->output_section_name, string("b_modes"), config->IsItBmodes);
		status = block->put_val<double>(config->output_section_name, string("theta_min"), config->thetamin);
		status = block->put_val<double>(config->output_section_name, string("theta_max"), config->thetamax);
		
		/*
		if(config->dolikelihood)
		{
			// make a resized data and covariance matrix to match the n_max
	        int nBands_in=int(config->BP_data.rows/nPairs);
	        //clog<<"nMaximum="<<nMaximum<<" n_max="<<n_max<<" nPairs"<<nPairs<<endl;
	        if(nBands>nBands_in)
	        {
	            clog<<"set nBands="<<nBands<<" is larger than the input nBands in the data vector:"<<nBands_in<<endl;
	            exit(1);
	        }
	        matrix BP_data(nBands*nPairs);
	        //En_data=config->En_data;
	        for(int i=0;i<nBands*nPairs;i++)
	        {
	            BP_data.load(i,config->BP_data.get((ceil(i/nBands)*nBands_in+i%nBands)));
	        }

	        matrix Cov_mat(nBands*nPairs,nBands*nPairs);
	        //Cov_mat=config->Cov_mat;
	        for(int i=0;i<nBands*nPairs;i++)
				for(int j=0;j<nBands*nPairs;j++)
					Cov_mat.load(i,j,config->Cov_mat.get((ceil(i/nBands)*nBands_in+i%nBands),(ceil(j/nBands)*nBands_in+j%nBands)));

			///Calculate likelihood here
			number ChiS=calChiS(BP_mat,BP_data,Cov_mat);
			number likelihood_val=-ChiS/2.0;
			string likename=config->output_section_name+"_like";
			//clog<<"likename is:"<<likename<<endl;
			status = block->put_val<number>(LIKELIHOODS_SECTION, likename, likelihood_val);
		}*/

	    return status;
	}
}// end of extern C


    
