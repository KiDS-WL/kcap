#ifndef KSICOVARIANCE_H
#define KSICOVARIANCE_H

#include "COSEBIs.h"

class KsiCovariance: public function_cosebis
{
public:

	KsiCovariance();
	KsiCovariance(number thetamin,number thetamax,int nKsiBins);
	~KsiCovariance();
	void setParameters(number thetamin1, number thetamax1,int nKsiBins1);
	void setParameters(matrix theta_plus, matrix theta_minus);
	///sets the noise parameters for calculating the covariance 
	void setNoise(number A1,vector<number> sigma_e_vec,vector<number> nBar_vec);
	///sets neff_vec
	void set_neff(vector<number> neff_vec);
	///sets noise value for each bin separately 
	void setNoise_vec(vector<number> noise_vec1);
	///sets noise to zero
	void SetNoNoise();
	//checks if input power soectrum is set. Returns false if not.
	bool checkPower();
	///set number of bin pairs and bins
	void setZbins(int nPairs1);
	///calculates power for a given redshift bin combination 
	void setPower(vector<number> log_ell_vec,vector<vector<number> > InputPower);
	///Return power
	number ReturnPower(number ell,int rPair);
	///sets 2PCFs from input KsiP and KsiM. Makes tables of them for interpolation. NOTE: the Ksi+- here need to
	///be from theory and smooth
	void setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM);
	///calculates 2PCFs for the given powerspectrum
	//void calKsi();
	///reads Npairs from an input ksi file from Athena
	void readNpairs(vector<string> FileName, int nColumns=8);
	///delta kronecker
	number delta(int i1, int j1);
	number integrant(number l);
	///calculates the elements of the covariance for a given theta1 and theta2, 
	///orderB=the order of the bessel functions used
	number value(number theta11,number theta21,int orderB11,int orderB21);
	vector<vector<number> > readKsi(string FileName,int nColumns);
	//vector<vector<number> > readfile(string FileName,int nColumns);
	vector<int> FindMinMaxTheta(vector<vector<number> > ksi_vec);
	void setIntegLimitTheta();
	///calculates the Ksi_pm covariance for the given parameters
	matrix CalCov(bool Logscale=true);
	matrix CalCovPlus(bool Logscale=true);
	///Calculate the mixed term of the covariance for xi_+, using a simple geometery and an area scaling
	matrix CalCovPlusMixed(bool Logscale=true);
	///calculates the number of redshift bin pair considered
	int calP(int nBins,int fbin,int sbin);
	// ///find the value of the covariance given a w2ww input file
	// number valueInputw2ww(int itheta1,int itheta2,char WhichKsi);
	// ///read in a w2ww file and find the range of theta that is needed
	// void readw2ww(vector<string> FileName_vec,int nColumns=8);
	// ///Calculate the mixed term of the covariance for xi_+, using an input w2ww calculated from data
	// matrix CalCovPlusMixed_Inputw2ww(bool Logscale=true);
	// ///uses the Cov in calCov to find a new covariance by interpolation, this doesn't work now
	// matrix CovInterpolate(int order,int NewNKsiBins,bool Logscale);

	
private:
	vector<function_cosebis> powerspectrum_vec;
  	vector<function_cosebis> Ksi_p_vec;
  	vector<function_cosebis> Ksi_m_vec;
  	vector <number> noise_vec,sigma_e_vec;
  	vector <matrix> Ksi_mat_vec;
 	vector <matrix> Tpm_mat_vec,pofz_mat_vec;
  	vector <matrix> Npair_mat_vec;
  	vector <number> neff_vec;
	vector<number> delta_vec;
	// vector<matrix> pofz_mat_vec;
	// vector<matrix> w2ww_mat_vec;
	// vector<KsiPlus> Ksi_p_vec;
 //  	vector<KsiMinus> Ksi_m_vec;
  	vector<number> Np;
  	vector<number> Theta_vec;
  	vector<int> integ_limit_theta;

	number A; //area

	number thetamin,thetamax,MinPowerCov,MaxPowerCov;

	int nBins,rp1,rp2,rp3,rp4,nPairs,redshiftPair;
	int nKsiBins;
	int orderB1,orderB2;
	number theta1,theta2,LogCoef;
//	int N_b; /// number of bins in w2ww
	int itheta_low;
	int N_b_range;

	number begin,end;
	number Noise,deltaNoise,noise,delta1noise,delta2noise,delta3noise,delta4noise;

	bool Logscale;
	bool NoNoise;
	bool MixedTerm;
};

#endif
