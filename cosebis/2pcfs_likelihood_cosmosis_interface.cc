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
	const string sectionName = SHEAR_XI_SECTION;
	const string shear_cl = SHEAR_CL_SECTION;

	///define a structure with everything that is needed to be read in setup and sent to execute
	// Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 

	struct pcfs_config {
		string input_section_name;
		string output_section_name;
		//int n_theta;
		number theta_min_plus;
		number theta_max_plus;
		number theta_min_minus;
		number theta_max_minus;
		vector<number> theta_plus;
		vector<number> theta_minus;
		//int input_option;//this can be either 0 for Cl or 1 for xi_pm
		string Input2PCFsFileName;
		string InputCovarianceFileName;
		//int noisy_input;//This is either 0: no noise or 1: noisy
		//COSEBIs *cosebis;
		matrix pcfs_data;
		matrix Cov_mat;
	} ;
	pcfs_config config;

	

	void *setup(cosmosis::DataBlock *options) 
	{
		//options reads the ini file
		//define config here and then read from options the relevant input quantities
// 		clog<<"I'm in status in COSEBIs_likelihood_cosmosis_inteface"<<endl;
  
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
	
		status=options->get_val<string>(sectionName, string("output_section_name"), string("2pcfs_results"), 
			config.output_section_name);
		status=options->get_val<number>(sectionName, string("theta_min_plus"), 1. , config.theta_min_plus);
		clog<<"got the value of theta_min_plus="<<config.theta_min_plus<<endl;
		status=options->get_val<number>(sectionName, string("theta_max_plus"),100., config.theta_max_plus);
		clog<<"got the value of theta_max_plus="<<config.theta_max_plus<<endl;
		status=options->get_val<number>(sectionName, string("theta_min_minus"), 1. , config.theta_min_minus);
		clog<<"got the value of theta_min_minus="<<config.theta_min_minus<<endl;
		status=options->get_val<number>(sectionName, string("theta_max_minus"),100., config.theta_max_minus);
		clog<<"got the value of theta_max_minus="<<config.theta_max_minus<<endl;
		// status=options->get_val<int>(sectionName, string("n_theta"),10, config.n_theta);
		// clog<<"got the value of n_theta="<<config.n_theta<<endl;

		status=options->get_val<string>(sectionName, string("input_2pcfs_filename"),"", config.Input2PCFsFileName);
		status=options->get_val<string>(sectionName, string("input_covariance_filename"),"", config.InputCovarianceFileName);
		matrix mat_in;
		///The input 2PCFS file should have 3 columns at least. 1: theta, 2: xip, 3: xim
		///The redshift pairs need to be in this order: 11, 12, 13, ... , nn
		matrix pcfs_in=mat_in.readFromASCII_marika((config.Input2PCFsFileName).c_str());
		matrix theta=pcfs_in.getColumn(0);
		matrix xip=pcfs_in.getColumn(1);
		matrix xim=pcfs_in.getColumn(2);
		///now make the input vector
		///Make cuts on PCFs based on the thetamin and thetamax given in the input.
		vector <number> theta_plus,theta_minus,xip_vec,xim_vec;
		for(int itheta=0; itheta<theta.rows; itheta++)
		{
			if((theta.get(itheta)>=config.theta_min_plus) && (theta.get(itheta)<=config.theta_max_plus))
			{
				theta_plus.push_back(theta.get(itheta)*arcmin);
				xip_vec.push_back(xip.get(itheta));
			}
			if((theta.get(itheta)>=config.theta_min_minus) && (theta.get(itheta)<=config.theta_max_minus))
			{
				theta_minus.push_back(theta.get(itheta)*arcmin);
				xim_vec.push_back(xim.get(itheta));
			}
		} 
		int ncols=theta_plus.size()+theta_minus.size();
		matrix pcfs_mat(1,ncols);
		for(int itheta=0; itheta<theta_plus.size(); itheta++)
		{
			pcfs_mat.load(0,itheta,xip_vec[itheta]);
		} 
		for(int itheta=0; itheta<theta_minus.size(); itheta++)
		{
			pcfs_mat.load(0,itheta+theta_plus.size(),xip_vec[itheta]);
		} 
		clog<<"theta_plus.size()="<<theta_plus.size()<<" theta_minus.size()="<<theta_minus.size()<<endl;
		config.pcfs_data=pcfs_mat;
		config.theta_plus=theta_plus;//change to radians
		config.theta_minus=theta_minus;//change to radians
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
		config.Cov_mat=mat_in.readFromASCII_marika((config.InputCovarianceFileName).c_str());
		status=options->get_val<string>(sectionName, string("input_section_name"), string("shear_xi"), 
			config.input_section_name);
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
// 		clog<<"I'm in execute in COSEBIs_likelihood_cosmosis_inteface"<<endl;

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

		vector<number> theta;
		status = block->get_val(config->input_section_name, string("theta"), theta);
		int ntheta=theta.size();
		
		if (status) 
		{
			clog<<"Could not load theta to 2PCFS liklihood"<<endl;
			return status;
		}

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
				string name_in_plus="xiplus_bin_"+toString(j_bin)+"_"+toString(i_bin);
				string name_in_minus="ximinus_bin_"+toString(j_bin)+"_"+toString(i_bin);
	    		status = block->get_val(config->input_section_name, name_in_plus, xip);
	    		status = block->get_val(config->input_section_name, name_in_minus, xim);
				if (status) 
				{
					clog<<"Could not load xip or xim bin "<<j_bin<<"_"<< i_bin<<" in 2PCFS likelihood"<<endl;
					return status;
				}
				
				Input_xip_vec_vec.push_back(xip);
				Input_xim_vec_vec.push_back(xim);
				///put Cl in a vector to send to COSEBIs.
				nPairs++;
			}
		}
		clog<<"setting 2PCFs"<<endl;
		
		
		///Make a vector of xip and xim based on the theta cuts and integrate over theta for the binning as well.
		matrix pcfs_mat=config->pcfs_data;
		int nTheta_minus=config->theta_minus.size()/nPairs;
		int nTheta_plus=config->theta_plus.size()/nPairs;
		vector<number> logtheta_minus(nTheta_minus),logtheta_plus(nTheta_plus);
		for (int itheta=0; itheta<nTheta_minus; itheta++)
			logtheta_minus[itheta]=log(config->theta_minus[itheta]);
		for (int itheta=0; itheta<nTheta_plus; itheta++)
			logtheta_plus[itheta]=log(config->theta_plus[itheta]);

		int counterP=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				function_cosebis pcfs_table;
				for(int r=0; r<Input_xip_vec_vec.size(); r++)
				{
					pcfs_table.loadWithValues(logtheta_plus,Input_xip_vec_vec[r],true);
					pcfs_table.extrapolationOff();
				}
				for(int itheta=0; itheta<nTheta_plus; itheta++)
				{
					pcfs_mat.load(0,counterP, pcfs_table.value(config->theta_plus[itheta]));
					counterP++;
				}
			}
		}


		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<=num_z_bin_B; j_bin++) 
			{
				function_cosebis pcfs_table;
				for(int r=0; r<Input_xim_vec_vec.size(); r++)
				{
					pcfs_table.loadWithValues(logtheta_minus,Input_xim_vec_vec[r],true);
					pcfs_table.extrapolationOff();
				}
				for(int itheta=0; itheta<nTheta_minus; itheta++)
				{
					pcfs_mat.load(0,counterP, pcfs_table.value(config->theta_minus[itheta]));
					counterP++;
				}
			}
		}


   //      // make a resized data and covariance matrix to match the n_max
   //      int nMaximum=int(config->En_data.rows/nPairs);
   //      //clog<<"nMaximum="<<nMaximum<<" n_max="<<n_max<<" nPairs"<<nPairs<<endl;
   //      if(n_max>nMaximum)
   //      {
   //          clog<<"set n_max="<<n_max<<" is larger than the folder's n_max:"<<nMaximum<<endl;
   //          exit(1);
   //      }
        
   //      matrix En_data(n_max*nPairs);
   //      for(int i=0;i<n_max*nPairs;i++)
   //          En_data.load(i,config->En_data.get((ceil(i/n_max)*nMaximum+i%n_max)));
		
   //      matrix Cov_mat(n_max*nPairs,n_max*nPairs);
   //      for(int i=0;i<n_max*nPairs;i++)
			// for(int j=0;j<n_max*nPairs;j++)
			// 	Cov_mat.load(i,j,config->Cov_mat.get((ceil(i/n_max)*nMaximum+i%n_max),(ceil(j/n_max)*nMaximum+j%n_max)));
           
//         clog<<"Cov_mat has "<<Cov_mat.rows<<" rows and "<<Cov_mat.columns<<" columns"<<endl;
//         clog<<"En_mat has "<<En_mat.rows<<" rows and "<<En_mat.columns<<" columns"<<endl;
//         clog<<"En_data has "<<En_data.rows<<" rows and "<<En_data.columns<<" columns"<<endl;
		///Calculate likelihood here
		matrix Delta2pcfs=config->pcfs_data-pcfs_mat;
		matrix iCov=config->Cov_mat.inverse();
		matrix chiS=Delta2pcfs.t()*iCov*Delta2pcfs;
// 	clog<<"chiS.size()="<<chiS.size()<<endl;
// 	clog<<"chiS.get(0)="<<chiS.get(0)<<endl;

		number ChiS=chiS.get(0);
		number likelihood_val=-ChiS/2.0;
// 		clog<<"likelihood="<<likelihood_val<<endl;
		status = block->put_val<number>(LIKELIHOODS_SECTION, "2pcfs_like", likelihood_val);
		return status;
	}
}// end of extern C



    
