#include "BandPower.h"


// constructor for the weight function calculator of band power
BandPower::BandPower()
{
	thetamin=1.;
	thetamax=1.;
	Response_function_type="tophat";

	gSet=false;
	Analytic=false;
	real=false;
	nBands=1;
	nPairs=1;
	FolderName="./cosebis/BandPower/";
	gFileName="g";
	WFileName="W";
	
	LLOW=1.;
	LHIGH=1e5;
	NLBINS=200.;
	bessel_order=0;
}


BandPower::~BandPower(){}


BandPower::BandPower(number thetamin,number thetamax, string Response_function_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName,string gFileName, string WFileName, bool real)
{
	initialize(thetamin,thetamax,Response_function_type,
		l_min_vec,l_max_vec,
		bessel_order,noApodise,Delta_x, Analytic,
		LLOW,
		LHIGH,
		NLBINS,
		FolderName,gFileName,WFileName,real);
  
}

///////
void BandPower::initialize(number thetamin,number thetamax, string Response_function_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName,string gFileName, string WFileName,bool real)
{
	clog<<"in BandPower initialize ......"<<endl;
	set_table_values(LLOW,LHIGH,NLBINS);
	setBandPowerName(FolderName,gFileName,WFileName);
	setTheta(thetamin,thetamax,Delta_x);
	gSet=false;
	WSet=false;
	set_noApodise(noApodise);
	Set_Response_function(Response_function_type,l_min_vec,l_max_vec);
	//setAnalytic(Analytic);
	set_bessel_order(bessel_order);
	setReal(real);
	if(real)
	{
		set_g(Analytic);
	}
	else
	{
		set_W();
	}
	nPairs=1;
}

void BandPower:: setAnalytic(bool Analytic1)
{
	Analytic=Analytic1;
}

void BandPower::set_table_values(number LLOW1,number LHIGH1, int NLBINS1)
{
	LLOW=LLOW1;
	LHIGH=LHIGH1;
	NLBINS=NLBINS1;
}

void BandPower::setBandPowerName(string FolderName1,string gFileName1,string WFileName1)
{
	FolderName=FolderName1;
	gFileName=gFileName1;
	WFileName=WFileName1;
}



void BandPower::setTheta(number thetamin1,number thetamax1,number Delta_x1)
{
	//clog<<"setting thetamin,thetamax"<<endl;
	thetamin=thetamin1;
	thetamax=thetamax1;
	//clog<<"thetamin="<<thetamin<<" thetamax="<<thetamax<<endl;
	Delta_x=Delta_x1;
}



void BandPower::set_bessel_order(int bessel_order1)
{
	bessel_order=bessel_order1;
	if((bessel_order==0) || (bessel_order==2) || (bessel_order==4) )
	{
		//clog<<"you gave bessel_order="<<bessel_order<<endl;
	}
	else
	{
		clog<<"you gave bessel_order="<<bessel_order<<endl;
		clog<<"not a recognised w_type, please use either of these: 0 for J0, 2 for J2 or 4 for J4, exiting now ..."<<endl;
		exit(1);
	}
}




void BandPower::set_noApodise(bool noApodise1)
{
	noApodise=noApodise1;
}


void BandPower::Set_Response_function(string type, vector<number> l_min_vec1,
														vector<number> l_max_vec1)
{
	Response_function_type=type;
	l_min_vec=l_min_vec1;
	l_max_vec=l_max_vec1;
	if(l_min_vec.size()!=l_max_vec.size())
	{
		clog<<"!!!!WARNING!!! l_min_vec and l_max_vec have different sizes"<<endl;
	}
	nBands=l_max_vec.size();
	set_l_min_max(l_min_vec[0],l_max_vec[nBands-1]);
}


void BandPower::set_l_min_max(number l_min1,number l_max1)
{
	l_min=l_min1;
	l_max=l_max1;
}

void BandPower::set_g(bool Analytic)
{
	if(!gSet)
	{
		//clog<<"g not set setting now:"<<endl;
		//BandPower_g g();
		number thetamin_g=exp(log(thetamin)-Delta_x/2.0);
		number thetamax_g=exp(log(thetamax)+Delta_x/2.0);
		g_vec.clear();
		for(int i=0; i<nBands; i++)
			g_vec.push_back(BandPower_g());
		for(int i=0; i<nBands; i++)
		{
			g_vec[i].initialize(thetamin_g,thetamax_g, Response_function_type,l_min_vec,l_max_vec,LLOW,LHIGH,Analytic,FolderName,gFileName);
			//g_vec[i].setAnalytic(Analytic);
		}
		gSet=true;
	}
}

//need to set the bessel_order and change it as well
void BandPower::set_W(bool Analytic)
{
	if(!WSet)
	{
		//clog<<"W not set setting now:"<<endl;
		//BandPower_W W();
		W_vec.clear();
		clog<<"bessel_order="<<bessel_order<<endl;
		for(int i=0; i<nBands; i++)
			W_vec.push_back(BandPower_W());
		for(int i=0; i<nBands; i++)
		{
			W_vec[i].initialize(thetamin,thetamax,Response_function_type,
				l_min_vec,l_max_vec,
				bessel_order,noApodise,Delta_x, Analytic,
				LLOW,
				LHIGH,
				NLBINS,
				FolderName,gFileName,WFileName);
			//W_vec[i].setAnalytic(Analytic);
			W_vec[i].set(i,bessel_order);
		}
		WSet=true;
	}
}

void BandPower::setReal(bool real1)
{
	real=real1;
}
//down to here in initialize
///////

// integrands for W depend on if apodisation is on or not and if tophat is used
number BandPower::integrant(number ell)
{
	//clog<<"in integ"<<endl;
	number integ= 0.;
	if(real)
	{
		//clog<<"real"<<endl;
		number theta=ell;
		integ=theta*Apodise(theta)*g_vec[bin_index].value(theta)*power_corr_vec[redshift].value(theta);
		
	}
	else
	{
		//clog<<"Fourier"<<endl;
		integ=ell*W_vec[bin_index].value(ell)*power_corr_vec[redshift].value(ell);
		//clog<<ell<<" "<<integ<<endl;
		//return integ;
	}
	return integ;
}


// void BandPower::print_integrant(bool noApodise1,int bessel_order1, number ell1, int bin_index1)
// {
// 	noApodise=noApodise1;
// 	ell=ell1;
// 	int nInteg=10000;
// 	bin_index=bin_index1;
// 	matrix integrant_mat(2,nInteg);
// 	if(noApodise)
// 	{
// 		for(int i=0; i<nInteg; i++)
// 		{
// 			number ellp=exp(log(LLOW)+log(LHIGH/LLOW)/(nInteg-1.)*i);
// 			integrant_mat.load(0,i,ellp);
// 			integrant_mat.load(1,i,integrant(ellp));
// 		}
// 	}
// 	else
// 	{
// 		for(int i=0; i<nInteg; i++)
// 		{
// 			number theta=exp(log(thetamin)+log(thetamax/thetamin)/(nInteg-1.)*i);
// 			integrant_mat.load(0,i,theta);
// 			integrant_mat.load(1,i,integrant(theta));
// 		}
// 	}
	
// 	integrant_mat.printOut(string("integrant_Ap.ascii").c_str(),5);
// }


void BandPower::setZbins(int nPairs1)
{
	//clog<<"setting nPairs: "<<nPairs1<<endl;
	nPairs=nPairs1;
	//clog<<"nPairs set:"<<nPairs<<endl;
}

void BandPower::setInput(vector<number> log_x, vector<vector<number> > Input)
{
	//clog<<"in set input, nPairs=";
	//nPairs=Input.size();
	//clog<<nPairs<<endl;

	if(Input.size()== nPairs)
	{
		//clog<<"Don't panic, number of redshift bin pairs matches the Input Power vector"<<endl;
	}
	else
	{
		clog<<"Panic! number of redshift bin pairs DOES NOT match the Input Power vector, exiting now ..."<<endl;
		exit(1);
	}

	power_corr_vec.clear();
	//clog<<"in set input, power_corr_vec cleared"<<endl;
	
	for(int r=0; r<nPairs; r++)
		power_corr_vec.push_back(function_cosebis());
	//clog<<"in set input, nPairs="<<endl;
	for(int r=0; r<nPairs; r++)
	{
		power_corr_vec[r].loadWithValues(log_x,Input[r],true);
		power_corr_vec[r].extrapolationOff();
// 		string powerFileName=string("Power_")+toString(r);
//   		power_corr_vec[r].setName(powerFileName.c_str(),function_cosebis::NONAMECOUNTER);
//  		power_corr_vec[r].makeTable(0.15,100000.,200,true);
//  		power_corr_vec[r].saveTable();
	}
	//clog<<"set the input"<<endl;
}


void BandPower::setInput_single(vector<number> log_x, vector<number> Input)
{
	//clog<<"in set input, nPairs=";
	//nPairs=Input.size();
	//clog<<nPairs<<endl;
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	nPairs=1;
	power_corr_vec.clear();
	power_corr_vec.push_back(function_cosebis());
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	power_corr_vec[0].loadWithValues(log_x,Input,true);
	power_corr_vec[0].extrapolationOff();
}

void BandPower::setInput_single_withExtrapolation(vector<number> log_x, vector<number> Input)
{
	//clog<<"in set input, nPairs=";
	//nPairs=Input.size();
	//clog<<nPairs<<endl;
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	nPairs=1;
	power_corr_vec.clear();
	power_corr_vec.push_back(function_cosebis());
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	power_corr_vec[0].loadWithValues(log_x,Input,true);
	power_corr_vec[0].extrapolationOn();
}

number BandPower::ReturnPower(number ell,int rPair)
{
	return power_corr_vec[rPair].value(ell);
}

number BandPower::ReturnW(number ell,int bin_index)
{
	return W_vec[bin_index].value(ell);
}

number BandPower::value(int bin_index1)
{
	number result=0.;
	bin_index=bin_index1;
	//clog<<"in value:"<<integ_limits_vec[bin_index].size()<<endl;
	for(unsigned int i=0; (i+1)<integ_limits_vec[bin_index].size(); i++)
	{
		//clog<<"in value:"<<integ_limits_vec[bin_index][i]<<" "<<integ_limits_vec[bin_index][i+1];
		number res=gaussianIntegrate_gsl(*this,integ_limits_vec[bin_index][i],integ_limits_vec[bin_index][i+1],20);
		//clog<<", res="<<res<<"     ";
		result+=res;
	}
	//clog<<"result="<<result<<endl;
	number Ni=valueNi(bin_index);
	//clog<<"Ni="<<Ni<<endl;
	if(real)
		return result*2.*pi/Ni;
	else
		return result/Ni;
}

//Note for now it only works for tophat response functions
number BandPower::valueNi(int bin_index)
{
	if(Response_function_type=="tophat")
	{
		return log(l_max_vec[bin_index]/l_min_vec[bin_index]);
	}
	else
	{
		//Haven't done this one
		return 1.;
	}
}

matrix BandPower::calBP()
{
	//clog<<"in calBP:"<<nBands*nPairs<<endl;
	matrix BP(nBands*nPairs);
	//bessel_order=bessel_order1;
	if(real)
	{
		set_g(Analytic);
		determine_integration_limits_real();
	}
	else
	{
		set_W(Analytic);
		determine_integration_limits_Fourier();
	}
	//clog<<"in calBP 2"<<endl;
	int band_index=0;
	for(int r=0; r<nPairs; r++)
	{
		redshift=r;
		//clog<<"redshiftpair="<<redshiftPair<<endl;
		for(int b=nBands*r,band_index=0 ;b<nBands*(r+1) ;b++,band_index++)
		{
			//clog<<"b="<<b<<endl;
			BP.load(b,value(band_index));
		}
	}
	//clog<<"calculated BP"<<endl;
	return BP;
}


matrix BandPower::calBP(vector<int> index)
{
	//clog<<"in calBP:"<<nBands*nPairs<<endl;
	matrix BP(index.size());
	//bessel_order=bessel_order1;
	if(real)
	{
		set_g(Analytic);
		determine_integration_limits_real();
	}
	else
	{
		set_W(Analytic);
		determine_integration_limits_Fourier();
	}

	redshift=0;
	for(int b=0; b<index.size(); b++)
	{
		//clog<<"b="<<b<<endl;
		BP.load(b,value(index[b]-1));
	}
	return BP;
}

void BandPower::determine_integration_limits_real()
{
	const int Nbins = 10000;
	integ_limits_vec.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	// number tmin=thetamin+exp(Delta_x);
	// number tmax=thetamax-exp(Delta_x);
	//clog<<"/////////////////////////////////////////////////////////////////"<<endl;
	//clog<<"tmin="<<tmin/arcmin<<endl;
	//clog<<"tmax="<<tmax/arcmin<<endl;
	//clog<<"/////////////////////////////////////////////////////////////////"<<endl;
	//exit(1);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(thetamin)+log(thetamax/thetamin)/(Nbins-1.)*i));
	}
	for(int b=0; b<nBands; b++)
	{
		for(int i=0;i<Nbins;i++)
		{
			table.load(1,i,g_vec[b].value(table.get(0,i)));
		}
		vector<number> integ_limits;
		integ_limits.push_back(thetamin);
		//integ_limits.push_back(tmin);
		for(int i=1;i<Nbins-1;i++)
		{
			if ((table.get(1,i-1)<=table.get(1,i) && table.get(1,i+1)<=table.get(1,i))
				|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
			{
				integ_limits.push_back(table.get(0,i));
			}
		}
		//integ_limits.push_back(tmax);
		integ_limits.push_back(thetamax);
		integ_limits_vec.push_back(integ_limits);
	}
}


void BandPower::determine_integration_limits_Fourier()
{
	const int Nbins = 10000;
	integ_limits_vec.clear();

	vector<number> integ_y;
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(LLOW)+log(LHIGH/LLOW)/(Nbins-1.)*i));
	}
	for(int b=0; b<nBands; b++)
	{
		for(int i=0;i<Nbins;i++)
		{
			table.load(1,i,W_vec[b].value(table.get(0,i)));
		}
		vector<number> integ_limits;
		integ_limits.push_back(LLOW);
		//integ_y.push_back(G_mu(ell,l_min_vec[bin_index],bessel_order));
		for(int i=1;i<Nbins-1;i++)
		{
			if ((table.get(1,i-1)<=table.get(1,i) && table.get(1,i+1)<=table.get(1,i))
				|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
			{
				integ_limits.push_back(table.get(0,i));
				//integ_y.push_back(table.get(1,i));
			}
		}
		integ_limits.push_back(LHIGH);
		integ_limits_vec.push_back(integ_limits);
	}
}


//T(theta) apodises the band powers to reduce oscillations around the edges
//Here x=ln(\theta) and Delta_x is the log width of the apodisation.
//x_l and x_u never go beyond the thetamin and thetamax range given by the user.
number BandPower::Apodise(number theta)
{
	if(Delta_x<=0)
		return 1.;

	number x= log(theta);
	number x_l=log(thetamin);
	number x_u=log(thetamax);
	number l1_bound=x_l-Delta_x/2.;
	number l2_bound=x_l+Delta_x/2.;
	number u1_bound=x_u-Delta_x/2.;
	number u2_bound=x_u+Delta_x/2.;
	number result=0.;
	//clog<<"x="<<x<<" l1_bound="<<l1_bound<<" l2_bound="<<l2_bound<<" u1_bound="<<u1_bound<<" u2_bound="<<u2_bound<<endl;
	if((l1_bound<=x) && (x<l2_bound))
	{
		//clog<<"in apodise: x_l="<<x_l<<" x="<<x<<endl;
		result= pow(cos(pi/2.*((x-(x_l+Delta_x/2.))/Delta_x)),2);
	}
	else if((l2_bound<=x) && (x<u1_bound))
		result= 1.;
	else if( (u1_bound<=x ) && (x<u2_bound) )
	{
		//clog<<"in apodise: x_u="<<x_u<<" x="<<x<<endl;
		result= pow(cos(pi/2.*((x-(x_u-Delta_x/2.))/Delta_x)),2);
	}
	else
		result= 0.;

	return result;
}
