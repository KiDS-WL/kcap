///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "cosmosis/datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "cosmosis/datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear two point correlation, 2PCFs, functions to E/B - COSEBIs
*/
#include "COSEBIs.h"

//namespace COSEBIS_ {
extern "C" {
	const string sectionName = "cosebis";
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
		int IsItNoisy;//This is either 0: no noise or 1: noisy
		string Input2PCFsFileName;
		string Input2PCFsCorrFileName;
		int input_2pcfs_filename_start;
		int input_2pcfs_filename_end;
		COSEBIs *cosebis;
		int nColumn;
		bool corr;
		int nZBins;
		bool range;
		int CalculateCov;
		string suffix;
	} ;
	COSEBIs_config config;
	

  	void *setup(cosmosis::DataBlock *options) 
  	{
  		//options reads the ini file
  		//define config here and then read from options the relevant input quantities
  		// clog<<"?????????????????????????????????????????????"<<endl;
  		// clog<<"I'm in setup in COSEBIs_2PCFs_cosmosis_inteface"<<endl;
  		// clog<<"????????????????????????????????????????????"<<endl;
  
  		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    	const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
    	
  		status=options->get_val<string>(sectionName, string("output_section_name"), string("COSEBIs_results"), 
  			config.output_section_name);
  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config.theta_min);
  // 		if (status) 
		// {
		// 	clog<<"Could not load theta_min to COSEBIs"<<endl;
		// 	clog<<"setting it to the default value of"<<1.<<endl;
		// 	//exit(1);
		// }
		// else
		// {
		clog<<"got the value of theta_min="<<config.theta_min<<endl;
		// }
  		//if (status) return failure;
  		//clog<<"reading the value of theta_max"<<endl;
	    status=options->get_val<number>(sectionName, string("theta_max"),100., config.theta_max);
	 //    if (status) 
		// {
		// 	clog<<"Could not load theta_max to COSEBIs"<<endl;
		// 	//exit(1);
		// }
		// else
		// {
		clog<<"got the value of theta_max="<<config.theta_max<<endl;
		//}
	    status=options->get_val<int>(sectionName, string("n_max"),10, config.n_max);
	 //    if (status) 
		// 	clog<<"Could not load n_max to COSEBIs"<<endl;
		// else
		clog<<"got the value of n_max="<<config.n_max<<endl;
		//the default is no noise=0
		status=options->get_val<int>(sectionName, string("noisy_input"), 0 , config.IsItNoisy);
		status=options->get_val<int>(sectionName, string("CalculateCov"), 0 , config.CalculateCov);
		status=options->get_val<string>(sectionName, string("CovNameSuffix"), "" , config.suffix);
		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
	    cosebis->initialize(config.n_max,config.theta_min,config.theta_max,1);
	    //npair set to one for now, will be set seperately in execute to the correct value
		clog<<"input is set to 2PCFs"<<endl;
		cosebis->setTs(config.n_max);
		status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_xi"), config.input_section_name);
		if (status) 
			clog<<"Could not load input section name"<<endl;
		else
			clog<<"got input section name: "<<config.input_section_name<<endl;

		status=options->get_val<int>(sectionName, string("n_zbins"),0, config.nZBins);
		clog<<"number of redshift bins is:"<<config.nZBins<<endl;
		if(config.IsItNoisy)
		{
			clog<<"Input is a noisy 2pcfs file, now reading file from disk"<<endl;
			status=options->get_val<string>(sectionName, string("input_2pcfs_filename"), config.Input2PCFsFileName);
			if (status) 
			{
				clog<<"Could not read the input 2PCFs filename, exiting now"<<endl;
				exit(1);
			}
			else
				clog<<"input 2PCFs filename is:"<<config.Input2PCFsFileName<<endl;
			status=options->get_val<string>(sectionName, string("input_2pcfs_correction_filename"), config.Input2PCFsCorrFileName);
			if (status) 
			{
				clog<<"Could not read the input 2PCFs correction filename, setting correction to none"<<endl;
				config.corr=false;
			}
			else
			{
				clog<<"input 2PCFs correction filename is:"<<config.Input2PCFsCorrFileName<<endl;
				config.corr=true;
			}

			status=options->get_val<int>(sectionName, string("input_2pcfs_filename_start"), config.input_2pcfs_filename_start);
			if (status) 
			{
				clog<<"no start range for the input files"<<endl;
				config.range=false;
			}
			else
			{
				clog<<"input 2PCFs filename starts at:"<<config.input_2pcfs_filename_start<<endl;
				status=options->get_val<int>(sectionName, string("input_2pcfs_filename_end"), config.input_2pcfs_filename_end);
				if (status) 
				{
					clog<<"no end range for the input files"<<endl;
					config.range=false;
				}
				else
				{
					clog<<"input 2PCFs filename ends at:"<<config.input_2pcfs_filename_end<<endl;
					config.range=true;
				}
			}

			status=options->get_val<int>(sectionName, string("number_of_columns_in_input_files"), 7, config.nColumn);
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

	    clog<<"/////////////////\\\\\\\\\\\\\\\\\\\\\\"<<endl;
	    clog<<"I'm in execute in COSEBIs_2PCFs_cosmosis_inteface"<<endl;
	    clog<<"/////////////////\\\\\\\\\\\\\\\\\\\\\\"<<endl;

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
		
		///do things for a noisy input. The integration routine is: trap
		if(config->IsItNoisy)
		{
			///a range of files is given as input
			for(int isim=start; isim<end; isim++)
			{
				bool ShouldIwriteTofile=false;
				clog<<"isim="<<isim<<endl;
				// this is the input file name for the noisy case: config->Input2PCFsFileName;
				///I need to know how many columns this file has and which one is theta, xi_+ and xi_-
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
									+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin);
							else
								xiPM_name=config->Input2PCFsFileName
									+string("_nBins_")+toString(config->nZBins)
									+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin);
							clog<<"input xiPM_name is: "<<xiPM_name<<endl;
							string xiPMCorr_name="";
							if(config->corr)
							{
								if(config->range)
									xiPMCorr_name=config->Input2PCFsCorrFileName+toString(isim)
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin)+"_Bin"+toString(j_bin);
								else
									xiPMCorr_name=config->Input2PCFsCorrFileName+toString(isim)
										+string("_Bin")+toString(i_bin)+string("_Bin")+toString(j_bin);
								FileCorrName_vec.push_back(xiPM_name);
							}
							FileName_vec.push_back(xiPM_name);
						}
					}
				}
				matrix EnFromFile;
				if(config->corr)
				{
					if(CheckFilesExist(FileName_vec) & CheckFilesExist(FileCorrName_vec))
					{
						EnFromFile=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec,FileCorrName_vec,config->nColumn);//do correction
						ShouldIwriteTofile=true;
					}
					else
					{
						clog<<"one or more files don't exist for isim="<<isim<<" moving on to the next"<<endl;
					}
				}
				else
				{
					if(CheckFilesExist(FileName_vec))
					{
						EnFromFile=config->cosebis->calEn2PCFsFromInputKsi(FileName_vec,config->nColumn);// nocorrection
						ShouldIwriteTofile=true;
					}
					else
					{
						clog<<"one or more files don't exist for isim="<<isim<<" moving on to the next"<<endl;
					}
				}
				if(ShouldIwriteTofile)
				{
					string FolderName="cosmosis-standard-library/cosebis/SLICS/Sims/";
					name_En=FolderName+"En"+config->suffix+"_LOS"+toString(isim)
						+string("_nBins_")+toString(config->nZBins)
						+"_"+thetaRange;
					name_Bn=FolderName+"Bn"+config->suffix+"_LOS"+toString(isim)
						+string("_nBins_")+toString(config->nZBins)
						+"_"+thetaRange;
					
					clog<<"write to file"<<endl;
					En_mat=EnFromFile.getColumn(1);
					Bn_mat=EnFromFile.getColumn(4);
					En_mat.printOut(name_En.c_str(),10);
					Bn_mat.printOut(name_Bn.c_str(),10);

					if(config->CalculateCov)
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
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange;
									name_Bn=string("cosebis_BnFrom2PCFs_LOS")+toString(isim)
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange;
								}
								else
								{
									name_En=string("cosebis_EnFrom2PCFs")
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange;
									name_Bn=string("cosebis_BnFrom2PCFs")
										+string("_nBins_")+toString(config->nZBins)
										+string("_Bin")+toString(i_bin+1)+string("_Bin")+toString(j_bin+1)+"_"+thetaRange;
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
			if(config->CalculateCov)
			{
				string FolderName="cosmosis-standard-library/cosebis/SLICS/Sims/";
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
				name_En=string("En_mean")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange;
				name_Bn=string("Bn_mean")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange;
				//status = block->put_val<vector<double> >(config->output_section_name, name_En, En_mean);
				//status = block->put_val<vector<double> >(config->output_section_name, name_Bn, Bn_mean);
				En_mean.printOut((FolderName+"/En_mean"+config->suffix+"_nBins_"+toString(config->nZBins)
					+"_nSims"+toString(nSims)+"_"+thetaRange).c_str(),10);
				Bn_mean.printOut((FolderName+"/Bn_mean"+config->suffix+"_nBins_"+toString(config->nZBins)
					+"_nSims"+toString(nSims)+"_"+thetaRange).c_str(),10);
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
				name_En=string("CovSimsEn")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange;
				name_Bn=string("CovSimsBn")+config->suffix+string("_nBins_")+toString(config->nZBins)+"_nSims"+toString(nSims)+"_"+thetaRange;
				//status = block->put_val<vector<double> >(config->output_section_name, name_En, CovSimsEn);
				//status = block->put_val<vector<double> >(config->output_section_name, name_Bn, CovSimsBn);
				CovSimsEn.printOut((FolderName+"/CovSimsEn"+config->suffix+"_nBins_"+toString(config->nZBins)
					+"_nSims"+toString(nSims)+"_"+thetaRange).c_str(),10);
				CovSimsBn.printOut((FolderName+"/CovSimsBn"+config->suffix+"_nBins_"+toString(config->nZBins)
					+"_nSims"+toString(nSims)+"_"+thetaRange).c_str(),10);
				clog<<"ChiSBn.size()="<<ChiSBn.size()<<endl;
				clog<<"ChiSBn.get(0,0)="<<ChiSBn.get(0,0)<<endl;
				clog<<"ChiSEn.get(0,0)="<<ChiSEn.get(0,0)<<endl;
			}
		}
		else
		{
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
			for (int i_bin=1; i_bin<=config->nZBins; i_bin++) 
			{
				for (int j_bin=i_bin; j_bin<=config->nZBins; j_bin++) 
				{
					vector<number> xiP_vec,xiM_vec;
					string name_in="xiplus_"+toString(j_bin)+"_"+toString(i_bin);
		    		status = block->get_val(config->input_section_name, name_in, xiP_vec);
					if (status) 
					{
						clog<<"Could not load bin "<<j_bin<<"_"<< i_bin<<" in XiP to COSEBIs"<<endl;
						return status;
					}
					name_in="ximinus_"+toString(j_bin)+"_"+toString(i_bin);
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




    
