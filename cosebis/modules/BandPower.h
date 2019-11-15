#ifndef  BANDPOWER_H
#define  BANDPOWER_H

#include "function_cosebis.h"
#include "Integrate.h"
#include <gsl/gsl_sf_bessel.h>
#include "BandPower_W.h"



//set the bessel order for W_vec and g_vec and have them open in the memory

//This class caluclates Band power either from input power spectra or \xi_\pm.
//It uses weight functions W^i_\mu(\ell) and g^i_\mu(\ell) in Fourier and real space respectively.
//The weight functions are defined as:
// W^i_\mu(\ell)=\int_{\theta_{min}}^{\theta_{max}} \d\theta \theta T(\theta) J_\mu(\ell\theta)g^i_\mu(\theta)
// i is the band power bin. 
// T(\theta) is the apodisation function designed to avoid the hard edges of the integral
// The general expresion for g^i_\mu is given by:
// g^i_\mu(\theta)=\int_0^\infty \d\ell \ell S^i(\ell)J_\mu(\ell\theta)
// where S^i(\ell) is the band power response function for bin i.
// g is calculated in BandPower_g and W is set is BandPower_W class.

//The band powers can be written as
//BP_\mu^i=1/N^i \int_0^\infty \d\ell \ell C(\ell) W_\mu^i (\ell)
// or
//BP_\mu^i=2\pi/N^i \int_{\theta_{min}}^{\theta_{max}} \d\theta \theta T(\theta) g^i(\theta) corr_\mu(\theta)
// with
// N^i=\int \d\ell/\ell S^i(\ell)
// This definition is desgined such that the band power traces \ell^2 C(\ell) at log centre of the bin i. 

// The four cases that we are interested are:
//BP_{clustering}=BP_0

//BP_{lensing E modes}= 1/(2 N^i) \int_0^\infty \d\ell \ell [(W_0+W_4)C_E(\ell)+(W_0-W_4)C_B(\ell)]
// = 1/2*{BP_0(C_E)+BP_4(C_B)+BP_0(C_B)-BP_4(C_B)}
// in real space
// = \pi/N^i \int_{\theta_{min}}^{\theta_{max}} \d\theta \theta T(\theta) [g^i_0(\theta) \xi_+(\theta)+ g^i_4(\theta) \xi_-(\theta)]
// = 1/2*{BP_0(\xi_+)+BP_4(\xi_-)}

//BP_{lensing B modes}= 1/(2 N^i) \int_0^\infty \d\ell \ell [(W_0-W_4)C_E(\ell)+(W_0+W_4)C_B(\ell)]
// = 1/2*{BP_0(C_E)-BP_4(C_B)+BP_0(C_B)+BP_4(C_B)}
// in real space
//  = \pi/N^i \int_{\theta_{min}}^{\theta_{max}} \d\theta \theta T(\theta) [g^i_0(\theta) \xi_+(\theta)- g^i_4(\theta) \xi_-(\theta)]
// = 1/2*{BP_0(\xi_+)-BP_4(\xi_-)}

// BP_{GGL}= BP_2

// There are three weight functions, denoted by bessel_order in the code:
// 1) bessel_order= 0 : J0
// 2) bessel_order= 4 : J4
// 3) bessel_order= 2 : J2

//if we set the response function S^i(\ell) to be a top hat between \ell_{low_i} and \ell_{high_i} then:
// g^i_\nu(\theta) simplifies depending on the order of \nu.
//

class BandPower : public function_cosebis
{
public:
	//sets some default values for everything.
	BandPower();
	//calls initialize
	BandPower(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName="./BandPower/",
		string gFileName="g",
		string WFileName="W",
		bool real=false);

	~BandPower();
	////////////////////////////////////////////////////////////////////////////////
	//initializes everything
	void initialize(number thetamin,number thetamax, string Response_function_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName="./BandPower/",
		string gFileName="g", 
		string WFileName="W",
		bool real=false);

	//sets lmin lmax and number of l_bins for making tables of the weight functions
	void set_table_values(number LLOW1,number LHIGH1, int NLBINS1);
	///sets the name of the folder and the start of the file names for g and W
	void setBandPowerName(string FolderName,string gFileName,string WFileName);
	//sets the type of weight to be calculated. Options are +, - and ggl
	void set_bessel_order(int bessel_order1);
	///sets thetamin thetamax and delta_x for apodisation. 
	///delta_x not used if noApodise set to True
	void setTheta(number thetamin1,number thetamax1, number Delta_x);
	void set_noApodise(bool noApodise);
	void setReal(bool real1);
	///set response function type and l_min, lmax and number of band powers:nBands
	void Set_Response_function(string type, vector<number> l_min_vec1,vector<number> l_max_vec1);
	///sets the l_min to the first l_min in l_min_vec and l_max to the last l_max in l_max_vec
	void set_l_min_max(number l_min1,number l_max1);
	/// if noApodise=true sets g by calling BandPower_g_pm
	void set_g(bool Analytic=true);
	void set_W(bool Analytic=true);
	///sets results to analytic if Analytic=true, if not top hat sets it back to false.
	void setAnalytic(bool Analytic);
	///////////////////////////////////////////////////////////////////////////
	void setZbins(int nPairs);
	//sets the input power spectrum
	void setInput(vector<number> log_ell,vector<vector<number> > InputPower);
	//sets the input power spectrum for a single bin
	void setInput_single(vector<number> log_ell,vector<number> InputPower);
	//this is the case where either there are no redshift bins or the input is given as one bin pair
	//void setInput(vector<number> log_ell,vector<number> InputPower);
	//returns the input power for the redshift bin and ell given as inputs
	number ReturnPower(number ell,int rPair);
	number ReturnW(number ell,int bin_index);
	//takes the bin index and calculates the Bandpower
	number value(int bin_index);
	///returns a matrix that has the band power values for the given bessel order and the input power or
	// correlation function if real is set to true
	matrix calBP();
	//only works for single bin, calculated BP for the list of indices given
	matrix calBP(vector<int> index);
	// vector<number> calBP_E();
	// vector<number> calBP_B();
	// vector<number> calBP_gg();
	// vector<number> calBP_ggl();
	//caluculates the normalisation for each band power
	number valueNi(int bin_index);
	//sets the power spectrum response function
	number S_response(int bin_index,number l);
	//writes response to file
	void writeResponse(int bin_index);
	/// integrant, depends on internal parameters: 
	/// (theta,bessel_order,bin_index) and real or Fourier
	number integrant(number x);
	void print_integrant(bool noApodise1,int bessel_order1, number ell1, int bin_index1);
	///determines integration limits for W
	void determine_integration_limits_Fourier();
	///determines integration limits for W without apodisation
	void determine_integration_limits_real();
	number Apodise(number theta);

private:
	/// internal parameters
	vector<BandPower_W> W_vec;
	vector<BandPower_g> g_vec;
	vector<number> l_min_vec,l_max_vec;
	vector<vector<number> > integ_limits_vec;
	vector<function_cosebis> power_corr_vec;

	number   thetamin;
	number   thetamax;
	number Delta_x;
	

	number LLOW;
	number LHIGH;
	int NLBINS;

	number l_min;
	number l_max;

	string FolderName,gFileName,WFileName;
	string Response_function_type;

	int bin_index,bessel_order,redshift;
	int nPairs;
	int nBands,nBins;

	bool noApodise,gSet,WSet, Analytic, real;


};

#endif
