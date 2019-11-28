///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear two point correlation, 2PCFs, 
functions to E/B - COSEBIs
*/
#include "COSEBIs.h"


extern "C" {
	const string sectionName = "cosebis";
	const string shear_cl = SHEAR_CL_SECTION;
	const number MinPowerCOSEBIs=0.1;
	const number MaxPowerCOSEBIs=1e6;
	const int PowerTableNumberCOSEBIs=200;

	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 

	typedef struct COSEBIs_config {
		string sectionName;
		string input_section_name_plus;
		string input_section_name_minus;
		string output_section_name;
		int n_max;
		number theta_min;
		number theta_max;
		int IsItNoisy;//This is either 0: no noise or 1: noisy
		string Input2PCFsFileName;
		string Input2PCFsCorrFileName;
		string Input2PCFsFileNameEnd;
		int input_2pcfs_filename_start;
		int input_2pcfs_filename_end;
		COSEBIs *cosebis;
		//int nColumn;
		bool corr;
		int nZBins;
		bool range;
		int CalculateCov;
		string suffix;
		bool calculate_c1_c2;
		string xipm_c1_c2_1;
		string xipm_c1_c2_2;
		number c1_1;
		number c2_1;
		number c1_2;
		number c2_2;
		string FolderName;
		vector<number> Correction_vec;
	} COSEBIs_config;
	
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
		
  		status=options->get_val<string>(sectionName, string("output_section_name"), string("cosebis"), 
  			config->output_section_name);
  		clog<<"output section name="<<config->output_section_name<<endl;

  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config->theta_min);
  		clog<<"theta_min="<<config->theta_min<<endl;

  		status=options->get_val<number>(sectionName, string("theta_max"),100., config->theta_max);
  		clog<<"theta_max="<<config->theta_max<<endl;

  		status=options->get_val<int>(sectionName, string("n_max"),10, config->n_max);
		clog<<"n_max="<<config->n_max<<endl;

  		status=options->get_val<string>(sectionName, string("OutputFolderName"), string("."), config->FolderName);
		clog<<"setting output folder name to "<<config->FolderName<<endl;

		status=options->get_val<int>(sectionName, string("CalculateCov"), 0 , config->CalculateCov);
		if(config->CalculateCov)
			clog<<"Will calculate covariance if possible"<<endl;

		status=options->get_val<string>(sectionName, string("NameSuffix"), "" , config->suffix);

		

		// read in xi_pm files with Athena format. These are files created using a constant c1 and c2 catalouge with the
		// same positions as the catalouge to be analysed. 
		// Simply swap the epsilon_1 and epsilon_2 in your input cats with c1 and c2 respectively. 
		// Do this for two sets of c1 and c2 so that the code can calculate sum (cos 4phi) and sum (sin 4phi) from the xi_mimus
		// component. This is used for the constant c-term modelling which affects xi_minus given 
		// a finite number of galaxies and  a finite field. 
		// The important part is that the first three columns are
		// theta xi_plus xi_minus
		// 
		string xipm_c1_c2_1,xipm_c1_c2_2;
		number c1_1, c1_2,c2_1, c2_2;
		status=options->get_val<string>(sectionName, string("xipm_c1_c2_1") , config->xipm_c1_c2_1);
		if(status)
		{
			clog<<"no input xipm_c1_c2_1 given"<<endl;
			config->calculate_c1_c2=false;
		}
		else
		{
			status=options->get_val<string>(sectionName, string("xipm_c1_c2_2") , config->xipm_c1_c2_2);
			if(status)
			{
				clog<<"no input xipm_c1_c2_1 given"<<endl;
				config->calculate_c1_c2=false;
			}
			else
			{
				status=options->get_val<number>(sectionName, string("c1_1") , config->c1_1);
				if(status)
				{
					clog<<"no input c1_1 given"<<endl;
					config->calculate_c1_c2=false;
				}
				else
				{
					status=options->get_val<number>(sectionName, string("c2_1") , config->c2_1);
					if(status)
					{
						clog<<"no input c2_1 given"<<endl;
						config->calculate_c1_c2=false;
					}
					else
					{
						status=options->get_val<number>(sectionName, string("c1_2") , config->c1_2);
						if(status)
						{
							clog<<"no input c1_2 given"<<endl;
							config->calculate_c1_c2=false;
						}
						else
						{
							status=options->get_val<number>(sectionName, string("c2_2") , config->c2_2);
							if(status)
							{
								clog<<"no input c2_2 given"<<endl;
								config->calculate_c1_c2=false;
							}
							else
							{
								clog<<"I've got the values for xipm_c1_c2_1:"<<config->xipm_c1_c2_1<<" and xipm_c1_c2_2:"<<config->xipm_c1_c2_2<<endl;
								clog<<"c1_1="<<config->c1_1<<" c2_1="<<config->c2_1<<" c1_2="<<config->c1_2<<" c2_2="<<config->c2_2<<endl;
								clog<<"going to calculate E_cos4phi and E_sin4phi given the input xi_- values"<<endl;
								config->calculate_c1_c2=true;
							}
						}
					}
				}
			}
		}

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

		//  initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
	    cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1,WnFolderName,TnFolderName,OutputTnFolderName);
	    //npair set to one for now, will be set seperately in execute to the correct value
		clog<<"input is set to 2PCFs"<<endl;
		cosebis->setTs(config->n_max);
		status=options->get_val<string>(sectionName, string("input_section_name_plus"), string("shear_xi_plus"), config->input_section_name_plus);
		clog<<"input section name plus is: "<<config->input_section_name_plus<<endl;
		status=options->get_val<string>(sectionName, string("input_section_name_minus"), string("shear_xi_minus"), config->input_section_name_minus);
		clog<<"input section name minus is: "<<config->input_section_name_minus<<endl;


		status=options->get_val<int>(sectionName, string("n_zbins"),0, config->nZBins);
		clog<<"number of redshift bins is:"<<config->nZBins<<endl;

		status=options->get_val<int>(sectionName, string("noisy_input"), 0 , config->IsItNoisy);

		if(config->IsItNoisy)
		{
			clog<<"Assuming that the input 2PCFs are noisy. Going to use trapezoidal integration over them"<<endl;
			clog<<"Now reading file from disk"<<endl;
			status=options->get_val<string>(sectionName, string("input_2pcfs_filename"), config->Input2PCFsFileName);
			if (status) 
			{
				clog<<"Could not read the input 2PCFs filename, exiting now"<<endl;
				exit(1);
			}
			else
				clog<<"input 2PCFs filename is:"<<config->Input2PCFsFileName<<endl;

			status=options->get_val<string>(sectionName, string("input_2pcfs_suffix"), config->Input2PCFsFileNameEnd);
			if (status) 
			{
				clog<<"Could not read the input 2PCFs filename end"<<endl;
				config->Input2PCFsFileNameEnd="";
			}
			else
				clog<<"input 2PCFs filename end is:"<<config->Input2PCFsFileNameEnd<<endl;


			status=options->get_val<string>(sectionName, string("input_2pcfs_correction_filename"), config->Input2PCFsCorrFileName);
			if (status) 
			{
				clog<<"Could not read the input 2PCFs correction filename, setting correction to none"<<endl;
				config->corr=false;
			}
			else
			{
				clog<<"input 2PCFs correction filename is:"<<config->Input2PCFsCorrFileName<<endl;
				config->corr=true;
			}

			vector<number> m_vec;
			//NOTE: if m_vec size doesn't match the number of redshift bins correction will not be applied
			status=options->get_val<vector<number> >(sectionName, string("m_bias_values"), m_vec);
			if(status)
			{
				clog<<"no input m_vec given"<<endl;
			}
			else
			{
				clog<<"got the values of m_vec for "<<m_vec.size()<<" bins"<<endl;
				if(m_vec.size()!=config->nZBins)
				{
					clog<<"m_vec size doesn't match the number of redshift. Will not apply m_correction"<<endl;
				}
				else
				{
					for(int bin1=0; bin1<m_vec.size(); bin1++)
						clog<<"m_vec["<<bin1<<"]="<<m_vec[bin1]<<endl;

					for(int bin1=0; bin1<m_vec.size(); bin1++)
						for(int bin2=bin1; bin2<m_vec.size(); bin2++)
							config->Correction_vec.push_back((1.+m_vec[bin1])*(1.+m_vec[bin2]));
				}
			}

			status=options->get_val<int>(sectionName, string("input_2pcfs_filename_start"), config->input_2pcfs_filename_start);
			if (status) 
			{
				clog<<"no start range for the input files"<<endl;
				config->range=false;
			}
			else
			{
				clog<<"input 2PCFs filename starts at:"<<config->input_2pcfs_filename_start<<endl;
				status=options->get_val<int>(sectionName, string("input_2pcfs_filename_end"), config->input_2pcfs_filename_end);
				if (status) 
				{
					clog<<"no end range for the input files"<<endl;
					config->range=false;
				}
				else
				{
					clog<<"input 2PCFs filename ends at:"<<config->input_2pcfs_filename_end<<endl;
					config->range =true;
				}
			}
		}
		else
			clog<<"Assuming that the input 2PCFs are not noisy. Going to use Gauss-Legendre integration over them"<<endl;

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
	   /* int num_z_bin_A;
		int num_z_bin_B;

		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
		else
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			num_z_bin_B = num_z_bin_A;
		}*/
	    string thetaRange=toString(config->theta_min,2)+string("-")+toString(config->theta_max,2);
	    clog<<"thetaRange is: "<<thetaRange<<endl;
		matrix En_mat,Bn_mat;
		vector<matrix> En_mat_vec, Bn_mat_vec;
		vector<number> theta, logTheta;
		int n_max=config->n_max;
		vector<number> En_vec(n_max),Bn_vec(n_max);
		vector<int> n_vals(n_max);
		int start= 0;
		int end= 0;
		string name_En,name_Bn;
		if(config->range)
		{
			start=config->input_2pcfs_filename_start;
			end=config->input_2pcfs_filename_end+1;
		}
		else
		{
			start=0;
			end=1;
		}
		string xiPM_name="";
		string xiPMCorr_name="";

		//if calculate_c1_c2 
		if(config->calculate_c1_c2)
		{
			vector<string> FileName_vec;
			matrix En_c1_c2_1,En_c1_c2_2;
			for (int i_bin=1; i_bin<=config->nZBins; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=config->nZBins; j_bin++) 
				{
					string xiPM_name="";
					xiPM_name=config->xipm_c1_c2_1
							+string("_nBins_")+toString(config->nZBins)
							+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+string(".ascii");
					clog<<"input xipm_c1_c2_1 is: "<<xiPM_name<<endl;
					FileName_vec.push_back(xiPM_name);
				}
			}
			bool filesExits=true;
			if(CheckFilesExist(FileName_vec))
			{
				En_c1_c2_1=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec);
				filesExits=true;
			}
			else
			{
				clog<<"one or more files don't exist"<<endl;
				filesExits=false;
			}


			FileName_vec.clear();
			for (int i_bin=1; i_bin<=config->nZBins; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=config->nZBins; j_bin++) 
				{
					string xiPM_name="";
					xiPM_name=config->xipm_c1_c2_2
							+string("_nBins_")+toString(config->nZBins)
							+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+string(".ascii");
					clog<<"input xipm_c1_c2_2 is: "<<xiPM_name<<endl;
					FileName_vec.push_back(xiPM_name);
				}
			}
			if(CheckFilesExist(FileName_vec))
			{
				En_c1_c2_2=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec);// nocorrection
				filesExits &=true;
			}
			else
			{
				clog<<"one or more files don't exist"<<endl;
				filesExits=false;
			}

			//if all xipm_c1_c2_1 and xipm_c1_c2_2 files exist then makes the COSEBIs for c1, c2
			if(filesExits)
			{
				// Intminus= int dt t T_n(t) xi_m(t) is the third column of the output
				// therefore we need to devide it by two to get the contribution 
				// to COSEBIs correctly
				matrix En_c1_c2_1_minus=0.5*En_c1_c2_1.getColumn(3);
				matrix En_c1_c2_2_minus=0.5*En_c1_c2_2.getColumn(3);
				matrix En_sin4phi,En_cos4phi;

				number a1=config->c1_1*config->c1_1-config->c2_1*config->c2_1;
				number a2=2.*config->c1_1*config->c2_1;

				number b1=config->c1_2*config->c1_2-config->c2_2*config->c2_2;
				number b2=2.*config->c1_2*config->c2_2;

				//eqs are: 
				// a1*En_cos4phi+a2*En_sin4phi=En_c1_c2_1_minus
				// b1*En_cos4phi+b2*En_sin4phi=En_c1_c2_2_minus
				// 
				if(a1==0)
				{
					En_sin4phi=En_c1_c2_1_minus/a2;
					En_cos4phi=(En_c1_c2_2_minus-b2*En_sin4phi)/b1;
				}
				else if(b1==0)
				{
					En_sin4phi=En_c1_c2_2_minus/b2;
					En_cos4phi=(En_c1_c2_1_minus-a2*En_sin4phi)/a1;
				}
				else
				{
					En_sin4phi=(En_c1_c2_2_minus-b1/a1*En_c1_c2_1_minus)/(b2-b1/a1*a2);
					En_cos4phi=(En_c1_c2_1_minus-a2*En_sin4phi)/a1;
				}
				clog<<"got the values for sin and cos"<<endl;
				En_sin4phi.printOut((config->FolderName+"/En_sin4phi"
					+"_nBins_"+toString(config->nZBins)+"_"+thetaRange+string(".ascii")).c_str(),10);

				En_cos4phi.printOut((config->FolderName+"/En_cos4phi"
					+"_nBins_"+toString(config->nZBins)+"_"+thetaRange+string(".ascii")).c_str(),10);
			}
			else
			{
				clog<<"xi_c1_c2 files don't exits!!!!"<<endl;
			}
		} //end if calculate_c1_c2 
		
		///do things for a noisy input. The integration routine is: trap
		if(config->IsItNoisy)
		{
			///a range of files is given as input
			for(int isim=start; isim<end; isim++)
			{
				bool ShouldIwriteTofile=false;
				clog<<"isim="<<isim<<endl;
				// this is the input file name for the noisy case: config->Input2PCFsFileName;
				///I need to know which columns are theta, xi_+ and xi_-
				vector<string> FileCorrName_vec, FileName_vec;
				///no automatic bin names given as input
				if (config->nZBins<=0)
				{
					if(config->range)
						xiPM_name=config->Input2PCFsFileName+toString(isim);
					else //no range was given
						xiPM_name=config->Input2PCFsFileName;
					clog<<"input xiPM_name is: "<<xiPM_name<<endl;
					if(config->corr)
					{
						if(config->range)
							xiPMCorr_name=config->Input2PCFsCorrFileName+toString(isim);
						else //no range was given
							xiPMCorr_name=config->Input2PCFsCorrFileName;
					}
					FileName_vec.push_back(xiPM_name);
				}
				else
				{
					//if nZBins>0
					for (int i_bin=1; i_bin<=config->nZBins; i_bin++) 
					{
						for (int j_bin=i_bin; j_bin<=config->nZBins; j_bin++) 
						{
							string xiPM_name="";
							if(config->range)
								xiPM_name=config->Input2PCFsFileName+toString(isim)
									+string("_nBins_")+toString(config->nZBins)
									+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+config->Input2PCFsFileNameEnd;
							else
								xiPM_name=config->Input2PCFsFileName
									+string("_nBins_")+toString(config->nZBins)
									+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+config->Input2PCFsFileNameEnd;
							clog<<"input xiPM_name is: "<<xiPM_name<<endl;
							FileName_vec.push_back(xiPM_name);

							string xiPMCorr_name="";
							if(config->corr)
							{
								if(config->range)
									xiPMCorr_name=config->Input2PCFsCorrFileName+toString(isim)
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+config->Input2PCFsFileNameEnd;
								else
									xiPMCorr_name=config->Input2PCFsCorrFileName+toString(isim)
										+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin)+config->Input2PCFsFileNameEnd;
								FileCorrName_vec.push_back(xiPMCorr_name);
							}
						}
					}
				}// end of nZbins>0

				matrix EnFromFile;
				//if correction is set to true
				if(config->corr)
				{
					if(CheckFilesExist(FileName_vec) && CheckFilesExist(FileCorrName_vec))
					{
						EnFromFile=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec,FileCorrName_vec);//do correction
						ShouldIwriteTofile=true;
					}
					
					else
					{
						clog<<"one or more files don't exist for isim="<<isim<<" moving on to the next"<<endl;
					}
				}
				else if(CheckFilesExist(FileName_vec) && (config->Correction_vec.size()==FileName_vec.size()))
				{
						clog<<"m_bias values read from ini are used to correct the xi_pm"<<endl;
						EnFromFile=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec,config->Correction_vec);//do correction
						ShouldIwriteTofile=true;
				}
				else if(CheckFilesExist(FileName_vec))
				{
					EnFromFile=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec);// nocorrection
					ShouldIwriteTofile=true;
				}
				else
				{
					clog<<"one or more files don't exist for isim="<<isim<<" moving on to the next"<<endl;
				}

				if(ShouldIwriteTofile)
				{
					if(config->range)
					{
						name_En=config->FolderName+"En"+config->suffix+"_LOS"+toString(isim)
							+string("_nBins_")+toString(config->nZBins)
							+"_"+thetaRange;
						name_Bn=config->FolderName+"Bn"+config->suffix+"_LOS"+toString(isim)
							+string("_nBins_")+toString(config->nZBins)
							+"_"+thetaRange;
					}
					else
					{
						name_En=config->FolderName+"En"+config->suffix
							+string("_nBins_")+toString(config->nZBins)
							+"_"+thetaRange+string(".ascii");
						name_Bn=config->FolderName+"Bn"+config->suffix
							+string("_nBins_")+toString(config->nZBins)
							+"_"+thetaRange+string(".ascii");
					}
					
					clog<<"write to file"<<endl;
					En_mat=EnFromFile.getColumn(1);
					Bn_mat=EnFromFile.getColumn(4);
					En_mat.printOut(name_En.c_str(),10);
					Bn_mat.printOut(name_Bn.c_str(),10);

					if(config->CalculateCov && config->range)
					{
						En_mat_vec.push_back(En_mat);
						Bn_mat_vec.push_back(Bn_mat);
					}
					//now write to file
					int p1=0;
					if (config->nZBins<=0)
					{
						for(int n=0;n<n_max;n++)
						{
							En_vec[n]=En_mat.get(n);
							Bn_vec[n]=Bn_mat.get(n);
						}
						if(config->range)
						{
							clog<<"range is true"<<endl;
							name_En=string("cosebis_EnFrom2PCFs_LOS")+toString(isim);
							name_Bn=string("cosebis_BnFrom2PCFs_LOS")+toString(isim);
						}
						else
						{
							name_En=string("cosebis_EnFrom2PCFs");
							name_Bn=string("cosebis_BnFrom2PCFs");
						}
						status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
						status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_vec);
					}
					else
					{

						for(int i_bin=0; i_bin<config->nZBins; i_bin++) 
						{
							for (int j_bin=i_bin; j_bin<config->nZBins; j_bin++) 
							{
								int m=0;
								//int p1=cosebis.calP(nBins,bin1,bin2);
								for(int n1=n_max*p1,m=0 ;n1<n_max*(p1+1) ;n1++,m++)
								{
									En_vec[m]=En_mat.get(n1);
									Bn_vec[m]=Bn_mat.get(n1);
								}
								if(config->range)
								{
									name_En=string("cosebis_EnFrom2PCFs_LOS")+toString(isim)
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange+string(".ascii");
									name_Bn=string("cosebis_BnFrom2PCFs_LOS")+toString(isim)
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange+string(".ascii");
								}
								else
								{
									name_En=string("cosebis_EnFrom2PCFs")
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange+string(".ascii");
									name_Bn=string("cosebis_BnFrom2PCFs")
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange+string(".ascii");
								}
								status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
								status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_vec);
								p1++;

							}
						}
					}
					clog<<"wrote COSEBIs to file for isim="<<isim<<endl;
				}
			}

			//calculate covariance
			if(config->CalculateCov && config->range)
			{
				
				int nSims=En_mat_vec.size();
				clog<<"now calculating the covariance matrix from "<<nSims<<" mocks"<<endl;
				matrix En_Sum=En_mat_vec[0];
				matrix Bn_Sum=Bn_mat_vec[0];

				for(int nSim=1; nSim<nSims; nSim++)
				{
					En_Sum+=En_mat_vec[nSim];
					Bn_Sum+=Bn_mat_vec[nSim];
				}
				
				matrix En_mean=En_Sum/number(nSims);
				matrix Bn_mean=Bn_Sum/number(nSims);
				name_En=string("En_mean")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange+string(".ascii");
				name_Bn=string("Bn_mean")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange+string(".ascii");
				//status = block->put_val<vector<double> >(config->output_section_name, name_En, En_mean);
				//status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_mean);
				En_mean.printOut((config->FolderName+name_En).c_str(),10);
				Bn_mean.printOut((config->FolderName+name_Bn).c_str(),10);
				matrix En2_Sum=En_mat_vec[0].t()*En_mat_vec[0].t()-En_mean.t()*En_mean.t();
				matrix Bn2_Sum=Bn_mat_vec[0].t()*Bn_mat_vec[0].t()-Bn_mean.t()*Bn_mean.t();
				
				for(int nSim=1; nSim<nSims; nSim++)
				{
					matrix En2=En_mat_vec[nSim].t()*En_mat_vec[nSim].t()-En_mean.t()*En_mean.t();
					matrix Bn2=Bn_mat_vec[nSim].t()*Bn_mat_vec[nSim].t()-Bn_mean.t()*Bn_mean.t();
					En2_Sum=En2_Sum+En2;
					Bn2_Sum=Bn2_Sum+Bn2;			
				}

				matrix CovSimsEn=(En2_Sum)/number(nSims-1.)/number(nSims);
				matrix CovSimsBn=(Bn2_Sum)/number(nSims-1.)/number(nSims);
				matrix varianceEn=CovSimsEn.diag();
				matrix varianceBn=CovSimsBn.diag();
				matrix stdEn=varianceEn.power(0.5);
				matrix stdBn=varianceBn.power(0.5);
				matrix iCovBn=CovSimsBn.inverse();
				matrix iCovEn=CovSimsEn.inverse();
				matrix ChiSBn=Bn_mean.t()*iCovBn*Bn_mean/number(Bn_mean.rows);
				matrix ChiSEn=En_mean.t()*iCovEn*En_mean/number(En_mean.rows);
				clog<<"En_mean.get(0,0)="<<En_mean.get(0)<<endl;
				clog<<"Bn_mean.get(0)="<<Bn_mean.get(0)<<endl;
				clog<<"CovSimsEn.get(0,0)="<<CovSimsEn.get(0,0)<<endl;
				clog<<"CovSimsBn.get(0,0)="<<CovSimsBn.get(0,0)<<endl;
				//clog<<"Saving coavriance to: "<<FolderName<<endl;
				name_En=string("CovSimsEn")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange+string(".ascii");
				name_Bn=string("CovSimsBn")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange+string(".ascii");
				//status = block->put_val<vector<double> >(config->output_section_name, name_En, CovSimsEn);
				//status = block->put_val<vector<double> >(config->output_section_name, name_Bn, CovSimsBn);
				CovSimsEn.printOut((config->FolderName+name_En).c_str(),10);
				CovSimsBn.printOut((config->FolderName+name_Bn).c_str(),10);
				clog<<"ChiSBn.size()="<<ChiSBn.size()<<endl;
				clog<<"ChiSBn.get(0,0)="<<ChiSBn.get(0,0)<<endl;
				clog<<"ChiSEn.get(0,0)="<<ChiSEn.get(0,0)<<endl;
			}
		}//end of noisy
		else // this needs to be updated because now xi_plus and minus are in different sections
		{
			status = block->get_val(config->input_section_name_plus, string("theta"), theta);
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
			for (int i_bin=1; i_bin<=config->nZBins; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=config->nZBins; j_bin++) 
				{
					vector<number> xiP_vec,xiM_vec;
					string name_in="bin_"+toString(j_bin)+"_"+toString(i_bin);
		    		status = block->get_val(config->input_section_name_plus, name_in, xiP_vec);
					if (status) 
					{
						clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in XiP to COSEBIs"<<endl;
						return status;
					}
					name_in="bin_"+toString(j_bin)+"_"+toString(i_bin);
		    		status = block->get_val(config->input_section_name_minus, name_in, xiM_vec);
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
			//fix this
			matrix En_all=config->cosebis->calEn2PCFs();
			En_mat=En_all.getColumn(1);
			Bn_mat=En_all.getColumn(4);

			//now write to file
			int p1=0;
			if (config->nZBins==0)
			{
				for(int n=0;n<n_max;n++)
				{
					En_vec[n]=En_mat.get(n);
					Bn_vec[n]=Bn_mat.get(n);
				}
				name_En=string("cosebis_EnFrom2PCFs");
				name_Bn=string("cosebis_BnFrom2PCFs");
				status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
				status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_vec);
			}
			else
			{
				for(int i_bin=0; i_bin<config->nZBins; i_bin++) 
				{
					for (int j_bin=i_bin; j_bin<config->nZBins; j_bin++) 
					{
						int m=0;
						//int p1=cosebis.calP(nBins,bin1,bin2);
						for(int n1=n_max*p1,m=0 ;n1<n_max*(p1+1) ;n1++,m++)
						{
							En_vec[m]=En_mat.get(n1);
							Bn_vec[m]=Bn_mat.get(n1);
						}
						name_En=string("cosebis_EnFrom2PCFs_")+string("_bin_")+toString(i_bin+1)+string("_")+toString(j_bin+1);
						name_Bn=string("cosebis_BnFrom2PCFs_")+string("_bin_")+toString(i_bin+1)+string("_")+toString(j_bin+1);
						status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
						status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_vec);
						p1++;

					}
				}
			}
		}

		for(int n=0;n<n_max;n++)
			n_vals[n]=n+1;
		status = block->put_val<vector<int> >(config->output_section_name, string("cosebis_n"), n_vals);
	    return status;
	}
}// end of extern C




    
