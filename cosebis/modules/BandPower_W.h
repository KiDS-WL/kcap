#ifndef  BANDPOWER_W_H
#define  BANDPOWER_W_H

#include "function_cosebis.h"
#include "Integrate.h"
#include <gsl/gsl_sf_bessel.h>
#include "BandPower_g.h"

//const int nTheta=100000;

//This class caluclates the weight functions for band powers and saves them into tables
//These are then called from the disk in the proceeding runs of the program. 
//The weight functions are defined as:
// W^i(ell)=\int_{\theta_{min}}^{\theta_{max}} \d\theta \theta T(\theta) J_\nu(\ell\theta)g^i_\nu(\theta)
// i is the band power bin. 
// T(\theta) is the apodisation function designed to avoid the hard edges of the integral
// The general expresion for g^i_\nu is given by:
// g^i_\nu(\theta)=\int_0^\infty \d\ell \ell S^i(\ell)J_\nu(\ell\theta)
// where S^i(\ell) is the band power response function for bin i such that given a 2D power spectrum C(\ell):
// C^i=1/N^i \int_0^\infty \d\ell \ell S^i(\ell) C(\ell)
// and
// N^i=\int \d\ell/\ell S^i(\ell)
// This definition is desgined such that the band power traces \ell^2 C(\ell) at log centre of the bin i. 

// There are three weight functions, denoted by bessel_order in the code:
// 1) bessel_order= 0 : J0
// 2) bessel_order= 4 : J4
// 3) bessel_order= 2 : J2

//if we set the response function S^i(\ell) to be a top hat between \ell_{low_i} and \ell_{high_i} then:
// g^i_\nu(\theta) simplifies depending on the order of \nu.
class BandPower_W : public function_cosebis
{
public:
	//sets some default values for everything.
	BandPower_W();
	//calls initialize
	BandPower_W(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName="./BandPower/",
		string gFileName="g",
		string WFileName="W");

	~BandPower_W();
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
		string WFileName="W");

	//sets lmin lmax and number of l_bins for making tables of the weight functions
	void set_table_values(number LLOW1,number LHIGH1, int NLBINS1);
	///sets the name of the folder and the start of the file names for g and W
	void setBandPower_WName(string FolderName,string gFileName,string WFileName);
	//sets the type of weight to be calculated. Options are +, - and ggl
	void set_bessel_order(int bessel_order1);
	///sets thetamin thetamax and delta_x for apodisation. 
	///delta_x not used if noApodise set to True
	void setTheta(number thetamin1,number thetamax1, number Delta_x);
	void set_noApodise(bool noApodise);
	///set response function type and l_min, lmax and number of band powers:nBands
	void Set_Response_function(string type, vector<number> l_min_vec1,vector<number> l_max_vec1);
	///sets the l_min to the first l_min in l_min_vec and l_max to the last l_max in l_max_vec
	void set_l_min_max(number l_min1,number l_max1);
	/// if noApodise=true sets g by calling BandPower_g_pm
	void set_g(bool Analytic);
	///sets results to analytic if Analytic=true, if not top hat sets it back to false.
	//void setAnalytic(bool Analytic);
	///////////////////////////////////////////////////////////////////////////

	//sets the power spectrum response function
	number S_response(int bin_index,number l);
	//writes response to file
	void writeResponse(int bin_index);
	// analytic solution for G given a top hat response and no apodisation
	number G_mu(number ell, number ellp,int mu);
	/// integrant for g_+-(theta), depends on internal parameters: 
	/// (theta,bessel_order,bin_index)
	number integrant(number x);
	void print_integrant(bool noApodise1,int bessel_order1, number ell1, int bin_index1);
	///set the bin_index and bessel order, then either load from file or calculate using get.
	void set(int bin_index1,int bessel_order1);
	///calculates g_pm for a given theta
	number get(number theta1);
	///returns bin_index value
	int show_bin_index();
	///determines integration limits for W
	void determine_integration_limits();
	///determines integration limits for W without apodisation
	void determine_integration_limits_W_noAp();
	///
	number W_noApodise(number ell);
	///returns the analytic value of g_pm if the resposnse function: S_i(ell) is a top hat.
	//number Theta_g_pm_TopHat(number theta);
	number Apodise(number theta);
	number Theta_g_tophat(number theta);

private:
	/// internal parameters
	BandPower_g g;
	vector<number> integ_limits,l_min_vec,l_max_vec;


	number   thetamin;
	number   thetamax;
	number   Delta_x; //This is?
	number   ell;
	int bin_index,bessel_order;

	number LLOW;
	number LHIGH;
	number lthresh; //switch integration method
	int NLBINS;

	number l_min;
	number l_max;

	string FolderName,gFileName,WFileName;
	string Response_function_type;

	bool noApodise,gSet, Analytic;
	int nBands,nBins,iband;

};

#endif
