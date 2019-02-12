///this class calculates En and Bn and their covariance matrix (Gaussian)
//It needs a convergence power spectrum or two point correlation functions as input
//for the covariance matrix it needs the noise value.

#ifndef COSEBIS_H
#define COSEBIS_H

//This one calculates logarithmic Wn functions which are used in making En from theory:
//En=int_0^inf dl l/(2pi) Wn(ell) P_E(ell)
//Bn=int_0^inf dl l/(2pi) Wn(ell) P_B(ell)
#include "WnLog.h"
//This one calculates T_+- which are used to calculate COSEBIs from 2PCFs
//En=int_{theta_min}^{theta_max} d theta theta/2[T_+(theta)xi_+(theta)+T_-(theta)xi_-(theta)]
//Bn=int_{theta_min}^{theta_max} d theta theta/2[T_+(theta)xi_+(theta)-T_-(theta)xi_-(theta)]
#include "TpnLogRadian.h"
//This makes matrices and does matrix stuff (most of the functions are disables because they used NR)
#include "matrix.h"

///minimum and maximum ell for the convergence power spectrum and the number of points in the table
///for the rest of the calculations a spline interpolation estimates the value of P(l) for the 
///desired l
///read these from cosmosis
//const number MinPowerCOSEBIs=0.1;
//const number MaxPowerCOSEBIs=1e6;
//const int PowerTableNumberCOSEBIs=200;

//
class COSEBIs : public function_cosebis
{
public:
	COSEBIs();
	~COSEBIs();
	COSEBIs(int nMax,number thetamin, number thetamax, int nPairs,string WnFolderName,string TnFolderName,string OutputTnFolderName);
	//void initialize(int nMax,number thetamin, number thetamax, int nPairs);
	void initialize(int nMaximum,number thetamin, number thetamax, int nPairs
		,string WnFolderName1="cosmosis-standard-library/cosebis/WnLog/"
		,string TnFolderName1="cosmosis-standard-library/cosebis/TLogsRootsAndNorms/"
		,string OutputTnFolderName="cosmosis-standard-library/cosebis/TpnLog");
	void setZbins(int nPairs1);
	///sets the COSEBIs parameters
	void setEparam(int nMax1,number thetamin1, number thetamax1);
	///initializes Wn_vec
	void setWns(int nMax);
	///initializes Tp_n_vec and Tm_n_vec
	void setTs(int nMaximum);
	///sets the noise parameters for calculating the covariance 
	void setNoise(number A1,number sigmaE1,number nBar1);
	///sets noise value for each bin separately 
	void setNoise_vec(vector<number> noise_vec1);
	///calculates power for a given redshift bin combination 
	void setPower(vector<number> log_ell_vec,vector<vector<number> > InputPower);
	///Return power
	number ReturnPower(number ell,int rPair);
	///sets 2PCFs from input KsiP and KsiM. Makes tables of them for interpolation. NOTE: the Ksi+- here need to
	///be from theory and smooth
	void setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM);
	///integrand for En or Cov depends if Cov_on is true or not
	number integrant(number l);
	///return the integrand for the power spectrum case
	matrix returnIntegrandForPowerCase(vector <number> ell_vec, int n, int redshift);
	///min max finder for En setWns uses this and determines the integration limits once and for all
	void determine_integration_limits_En();
	///min max finder for the plus part of the En from 2PCFs integral
	void determine_integration_limits_En_Tp();
	///min max finder for the minus part of the En from 2PCFs integral
	void determine_integration_limits_En_Tm();
	///the integrand used to find En from the Trapezoidal integration rutine. if noisyKsi is true
	///noise will be added to Ksi_+ and Ksi_- before integration
	number trapIntegrant(number theta);
	///finds the value E_n from 2PCFs
	matrix valueEn2PCFs(int n);
	///calculates En for all n from 2PCFs
	matrix calEn2PCFs();
	///calculates En for a given value of n
	number valueEn(int n);
	///calculate En frm P(ell)
	matrix calEn();
	//calculates the covariance for En Em, n=n1 m=m1
	number valueCov(int n1,int m1);
	///calculates Covariance matrix and returns it
	matrix calCov();
	///this calculates the B-mode covariance assuming that there is no P_B(l) and noise is Gaussian, it returns
	///1/(2piA)*\int dl l W_n W_m
	number valueBCov(int n1,int m1);
	matrix calBCov();
	///returns 1 if i==j, zero otherwise
	number delta(int i, int j);
	///Determines the min max values for integration for the Covariance
	void determine_integration_limitsCov();
	///takes n for En and the number of bins for the trapezoidal integral and returns the En values
	///using 2PCFs
	//matrix valueEn2PCFsTrap(int n,int nBinsTrap,number h);
	//matrix calEn2PCFsTrap(int nBinsTrap,bool noisyKsi1);
	///calculates the En from the input Ksi_mat with trapezoidal integration resturns a matrix:
	///0:(Int_p+Int_m)/2 1:Int_p 2:Int_m
	matrix valueEn2PCFsKsiInput(matrix& Ksi_mat,int m,number h);
	///reads a Ksi file with theta ksi+ and ksi- returns the values in a 2D vector 
	vector<vector<number> > readKsi(string FileName,int nColumns);
	///finds the minTheta and maxTheta which are the closest to thetamin and thetamax
	///return the index of the MinTheta=index[0] and MaxTheta=index[1]
	vector<int> FindMinMaxTheta(vector<vector<number> > ksi_vec);
	///evaluates Tpm for the given theta's in theta_mat and saves them in Tpm_mat_vec, does this only once
	void FindTnKsiInput(matrix theta_mat);
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName,int nColumns=8);	
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName, vector<string> corrFile,int nColumns=8);
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName, vector<number> Corr_vec,int nColumns=8);
	///sets the value of the index from the bin
	int calP(int nBins,int fbin,int sbin);
	///reads En from a file
	matrix readInputEn(string InputEnfileName);
	///reads covariance of COSEBIs from a file
	matrix readInputCovariance(string InputCovfileName);
	///returns Chi^2 from En_data, covarinace and Encal
	number CalChiS(matrix En_th,matrix En_data,matrix Cov_mat);

	///takes En as input and returns its numerical derivative, h is the step of integration
	/// size is the lenght of the En matrix (rows) and the derivative mathod depends on the 
	/// variable derivative= 2,4,6 ==> 3,5,7 point stencil methods
	//matrix Derivative(matrix& EnM,number h,int size);
	
	//void setParamValueDEn(int i,int par,number h);
	///calculates derivatives of En
	//matrix calDEn();
	/// takes DEn_max as input to find M
	//matrix calM(int i, int j, int n);
	/// takes DEn_max as input to find M, changes the order of variables. W0 comes last. only for flat cosmology with max 6 params
	//matrix calMW0Last(int i, int j, int nMax);
	///returns the parameters to their fiducial value depending on their order
	//void initializeParam();
	///sets the value of parameters of the second order derivative
	//void setParamValueDDEn(int i, int par1,int par2,number h,number k);
	///sets the value of parameters of the simple second order derivative
	//void setParamValueDDEnSimple(int i, int par1,int par2,number h,number k);
	///another way of getting the secnd order derivative of En
	//matrix SecondDerivative(matrix & En_mat,number h, number k);
	///calculates the fourpoint second order partial derivative of the input matrix
	//matrix fourPointSecondDerivative(matrix& En_mat,number h,number k);
	///calculates H a 3 index quantity which contains the second order partial derivatives of En
	//vector<vector< matrix> > calDDE();
	///calculates H a 3 index quantity which contains 
	///the simple second order partial derivatives of En
	//vector<vector< matrix> > calDDESimple();
	///calculates Q= H_{n \mu \nu} \Phi_\nu   
	//matrix calQ(matrix deltaPhi,int nMax);

	///sets parameter order to the default order
	//void setParamsOrderToDefault();
	///sets the order of parameters for derivatives
	//void setParamsOrder(vector<string> paramsOrder1);

	
private:
  string WnFolderName,TnFolderName,OutputTnFolderName;
  vector<WnLog> Wn_vec;
  vector<TpnLog> Tpn_vec; 
  vector<function_cosebis> powerspectrum_vec;
  vector<function_cosebis> Ksi_p_vec;
  vector<function_cosebis> Ksi_m_vec;
  //RandomGaussian random;  

  vector <number> integ_limits,integ_limitsTp,integ_limitsTm;
  vector <vector<number> > integ_limits_vec;
  vector <number> noise_vec;
  vector <matrix> Ksi_mat_vec;
  vector <matrix> Tpm_mat_vec,pofz_mat_vec;
  //vector <string> paramsOrder;
  //vector <vector<matrix> > DDEn_vecvec;
  //matrix DEn_mat;
  //matrix pofz_mat;
  //matrix random_matP;
  //matrix random_matM;
  
  number delta1noise,delta2noise,delta3noise,delta4noise;  
  number A,sigmaE,begin,end,nBar;
  number thetamin,thetamax,lthresh,LHIGH,NoiseKsi;
  
  int param,nPairs,powerswitch,nBins,nMaximum,derivative;
  int redshiftPair,rp1,rp2,rp3,rp4;
  int nW,mW,counter,nT,iRand;  

  ///default is false
  bool TpmNotDone,noisyKsi, BCov_on,DEn_calculated,
		OneParam,Cov_on,Real,realPlus,WnSet,TnSet,EnInteglimitSet;

  matrix En_data, Cov_mat;
};

#endif
