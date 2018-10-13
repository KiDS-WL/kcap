#ifndef  WNLOG_H
#define  WNLOG_H

#include "function_cosebis.h"
#include "Integrate.h"
#include "besselzero04.h"
#include <gsl/gsl_sf_bessel.h>

/// the lowest ell that is calculated
const number LLOW = 1.;
const number high =500.;
///the highest ell that is calculated
//	LHIGH=high*20./thetamax;

// The number of points in the Gauss-Legendre integration
const number accuracyG=40;

/// number of log-bins for Wn-table
const int    NLBINS  = 1000000;
// const int    NLBINS  = 10000;

class WnLog : public function_cosebis
{
public:

	WnLog();
	WnLog(number thetamin1,number thetamax1,int nMax,string WnFileName="WnLog/WnLog");
	~WnLog();

	/// integrant for Wn(l), depends on internal parameters (thetamin,thetamax,n)
	number integrant(number x);
	///sets thetamin thetamax and nMax to read the roots and the normalization from disk
	void setTheta(number thetamin1,number thetamax1,int nMax);
	///sets the Folder and the start of the WnLog file name. 
	void setWnLogName(string WnFileName);
	/// sets the internal parameters (n) to open Wn file or make it
	void set(int order);
	///returns n 
	int show_n();
	/// returns Wn(l) 
	number get(number x);
	///writes TnLog on a matrix and on disk
	void writeTnLog(int n1);
	///calculates the value of WnLog using an step by step integration between zeros of bessel function and TnLogs
	number valueStep(number l);
	void StepFinder();
	void writeStep(int n1);
	//number LHIGH;


private:

	/// TLog+n(x)
	number TnLog(number theta);
	/// internal parameters
	function_cosebis Wn_table;
	vector<number> j0minmax;
	vector< vector<number> > root;
	vector< vector<number> > rootTheta;
	vector <number> norm;
	vector <number> StepN;
	int      n;
	number   l;
	number   lthresh;
	number   thetamin;
	number   thetamax;
	number   B;
	number LHIGH;
	string WnFileName;
};

#endif
