#ifndef  BANDPOWER_G_H
#define  BANDPOWER_G_H

#include "function_cosebis.h"
#include "Integrate.h"
#include <gsl/gsl_sf_bessel.h>

const int nTheta=10000;


// This calculates g+/- using eq 7:
// g^i_\mu=\int_0^\infty dl l S^i(l)J_\mu(l\theta)
// S is the band power response function
//input theta should be in radians
class BandPower_g : public function_cosebis
{
public:

	BandPower_g();
	BandPower_g(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,number LLOW, number LHIGH
		,string FolderName="./BandPower/"
		,string FileName="g");
	void initialize(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,number LLOW, number LHIGH
		,string FolderName="./BandPower/"
		,string FileName="g");
	~BandPower_g();
	///sets the Folder and file names
	void setBandPower_gName(string FolderName1,string FileName1);
	///sets results to analytic if Analytic=true, if not top hat sets it back to false.
	void setAnalytic(bool Analytic);
	///sets thetamin thetamax
	void setTheta(number thetamin1,number thetamax1);
	void Set_Response_function(string type, vector<number> l_min_vec1,vector<number> l_max_vec1);
	void set_l_min_max(number LLOW1,number LHIGH1);
	number S_response(int bin_index,number l);
	void writeResponse(int bin_index);
	/// integrant for g_+-(theta), depends on internal parameters: 
	/// (theta,bessel_order,bin_index)
	number integrant(number x);
	///set the bin_index and bessel order, then either load from file or calculate using get.
	void set(int bin_index1,int bessel_order1);
	///returns bin_index value
	int show_bin_index();
	///determines integration limits for g_pm
	void determine_integration_limits();
	///calculates g_pm for a given theta
	number get(number theta1);
	///returns the analytic value of g_pm if the resposnse function: S_i(ell) is a top hat.
	number Theta_g_tophat(number theta);

private:
	/// internal parameters
	vector<number> integ_limits,l_min_vec,l_max_vec;

	number   thetamin;
	number   thetamax;
	number   theta;
	int bin_index,bessel_order;
	int nBands;

	number LLOW;
	number LHIGH;

	bool Analytic;

	string FolderName,FileName;
	string Response_function_type;

};

#endif
