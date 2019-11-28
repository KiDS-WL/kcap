#include "BandPower_W.h"


// constructor for the weight function calculator of band power
BandPower_W::BandPower_W()
{
	thetamin=1.;
	thetamax=1.;
	Response_function_type="tophat";

	gSet=false;
	Analytic=false;
	nBands=1;
	// lthresh=0.;
	FolderName="./cosebis/BandPower/";
	gFileName="g";
	WFileName="W";
	
	LLOW=1.;
	LHIGH=1e5;
	NLBINS=200.;
	bessel_order=0;

	//clog<<"const in BandPower_W.........."<<endl;
}


BandPower_W::~BandPower_W(){}


BandPower_W::BandPower_W(number thetamin,number thetamax, string Response_function_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName,string gFileName, string WFileName)
{
	initialize(thetamin,thetamax,Response_function_type,
		l_min_vec,l_max_vec,
		bessel_order,noApodise,Delta_x, Analytic,
		LLOW,
		LHIGH,
		NLBINS,
		FolderName,gFileName,WFileName);
  
}
///////
void BandPower_W::initialize(number thetamin,number thetamax, string Response_function_type,
		vector<number> l_min_vec,vector<number> l_max_vec,
		int bessel_order, bool noApodise, number Delta_x, bool Analytic,
		number LLOW,
		number LHIGH,
		int NLBINS,
		string FolderName,string gFileName, string WFileName)
{
	clog<<"in BandPower_W initialize ......"<<endl;
	set_table_values(LLOW,LHIGH,NLBINS);
	setBandPower_WName(FolderName,gFileName,WFileName);
	set_bessel_order(bessel_order);
	setTheta(thetamin,thetamax,Delta_x);
	gSet=false;
	set_noApodise(noApodise);
	Set_Response_function(Response_function_type,l_min_vec,l_max_vec);
	setAnalytic(Analytic);
	if(!noApodise)
		set_g();
}


void BandPower_W::set_table_values(number LLOW1,number LHIGH1, int NLBINS1)
{
	LLOW=LLOW1;
	LHIGH=LHIGH1;
	NLBINS=NLBINS1;
}

void BandPower_W::setBandPower_WName(string FolderName1,string gFileName1,string WFileName1)
{
	FolderName=FolderName1;
	gFileName=gFileName1;
	WFileName=WFileName1;
}

void BandPower_W::set_bessel_order(int bessel_order1)
{
	bessel_order=bessel_order1;
	if( (bessel_order==0) || (bessel_order==2) || (bessel_order==4) )
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


void BandPower_W::setTheta(number thetamin1,number thetamax1,number Delta_x1)
{
	//clog<<"setting thetamin,thetamax"<<endl;
	thetamin=thetamin1;
	thetamax=thetamax1;
	//clog<<"thetamin="<<thetamin<<" thetamax="<<thetamax<<endl;
	Delta_x=Delta_x1;
}


void BandPower_W::set_noApodise(bool noApodise1)
{
	noApodise=noApodise1;
}


void BandPower_W::Set_Response_function(string type, vector<number> l_min_vec1,
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


void BandPower_W::set_l_min_max(number l_min1,number l_max1)
{
	l_min=l_min1;
	l_max=l_max1;
}


void BandPower_W::set_g()
{
	if(!gSet)
	{
		//clog<<"g not set setting now:"<<endl;
		g.initialize(thetamin,thetamax, Response_function_type,l_min_vec,l_max_vec,LLOW,LHIGH,FolderName,gFileName);
		gSet=true;
	}
}

void BandPower_W::setAnalytic(bool Analytic)
{
	if(!gSet)
	{
		set_g();
	}
	g.setAnalytic(Analytic);
}

//down to here in initialize
///////

number BandPower_W::S_response(int bin_index,number l)
{
	if(Response_function_type=="tophat")
	{
		if((l_min_vec[bin_index]<=l) && (l<l_max_vec[bin_index]))
			return 1.;
		else
			return 0.;
	}
	else if(Response_function_type=="SomeOtherType_to_be_defined")
	{
		Analytic=false;
		return 10.;
	}
	else
	{
		Analytic=false;
		return 0;
	}
}

void BandPower_W::writeResponse(int bin_index)
{
	matrix S_mat(2,10000);
	for(int i=0; i<10000; i++)
	{
		number ell=exp(log(l_min)+log(l_max/l_min)/(10000-1)*i);
		S_mat.load(0,i,ell);
		S_mat.load(1,i,S_response(bin_index,ell));
	}
	S_mat.printOut((FolderName+string("/S_")+toString(bin_index)+string("_")
		+toString(nBins)+string("_")+toString(l_min_vec[0],2)
		+string("-")+toString(l_max_vec[nBins-1],2)+string(".ascii")).c_str(),8);
}

//analytic solution for G given a top hat response
//G_\mu(\ell,\ell')=\int_{\theta_{min}}^{\theta_{max}} \d\theta \theta J_\mu(\ell\theta)J_\mu(\ell'\theta)
number BandPower_W::G_mu(number ell, number ellp,int mu)
{
	number G_mu_max=thetamax/(ellp*ellp-ell*ell)*(
		 ellp*gsl_sf_bessel_Jn(mu+1,ellp*thetamax)*gsl_sf_bessel_Jn(mu,ell*thetamax)
		-ell*gsl_sf_bessel_Jn(mu,ellp*thetamax)*gsl_sf_bessel_Jn(mu+1,ell*thetamax));
	number G_mu_min=thetamin/(ellp*ellp-ell*ell)*(
		 ellp*gsl_sf_bessel_Jn(mu+1,ellp*thetamin)*gsl_sf_bessel_Jn(mu,ell*thetamin)
		-ell*gsl_sf_bessel_Jn(mu,ellp*thetamin)*gsl_sf_bessel_Jn(mu+1,ell*thetamin));
	return G_mu_max-G_mu_min;
}

// integrands for W depend on if apodisation is on or not and if tophat is used
number BandPower_W::integrant(number theta)
{
	number integ=0.;
	if(noApodise)
	{
		//This is a general solution when no apodisation is set
		number ellp=theta;
		integ=ellp*G_mu(ell,ellp,bessel_order)*S_response(bin_index,ellp);
	}
	else
	{
		// this is a general solution when apodisation is set
		integ=theta*Apodise(theta)*gsl_sf_bessel_Jn(bessel_order,ell*theta)*g.value(theta);
	}
	return integ;
}


void BandPower_W::print_integrant(bool noApodise1,int bessel_order1, number ell1, int bin_index1)
{
	noApodise=noApodise1;
	string ap_str;

	if(noApodise)
		ap_str="noAp";
	else
	{
		ap_str="Ap";
		ap_str+=toString(Delta_x,2);
	}

	ell=ell1;
	int nInteg=10000;
	bin_index=bin_index1;
	matrix integrant_mat(2,nInteg);
	g.set(bin_index,bessel_order);

	clog<<"thetamin="<<thetamin<<"thetamax="<<thetamax<<endl;

	if(noApodise)
	{
		for(int i=0; i<nInteg; i++)
		{
			number ellp=exp(log(LLOW)+log(LHIGH/LLOW)/(nInteg-1.)*i);
			integrant_mat.load(0,i,ellp);
			integrant_mat.load(1,i,integrant(ellp));
		}
	}
	else
	{
		for(int i=0; i<nInteg; i++)
		{
			number theta=exp(log(thetamin)+log(thetamax/thetamin)/(nInteg-1.)*i);
			integrant_mat.load(0,i,theta/arcmin);
			integrant_mat.load(1,i,integrant(theta));
		}
	}
	
	integrant_mat.printOut((string("integrant_")+ap_str+"_bin"+toString(bin_index)+string("_ell_")+toString(ell,2)+string("_")+toString(bessel_order)+string(".ascii")).c_str(),5);
}

///here is where the weights are set either by calculating them or by loading the table from disk
void BandPower_W::set(int bin_index1,int bessel_order1)
{
	//clog<<"set bin_index and bessel_order for BandPower_W"<<endl;
	bin_index=bin_index1;
	bessel_order=bessel_order1;
	// is there a table on disk?
	string pm_string=toString(bessel_order);

	string myname =FolderName+string("/")+WFileName+string("_")+pm_string+string("_")
			+toString(l_min_vec[bin_index],2)+string("-")+toString(l_max_vec[bin_index],2)+string("_")
		 	+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
	ifstream fhandler((myname+string(".table")).c_str());
	// NO?
	if (fhandler.fail())
	{
		if(!CheckFolderExist(FolderName))
		{
			//clog<<"making folder for BandPower_W:"<<FolderName<<endl;
			mkdir((FolderName).c_str(), 0777);
		}
		if(!CheckFolderExist(FolderName))
		{
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!! Can't make the Folder for Wight functions!!!!"<<endl;
			clog<<"!!!!!!!!!!! WILL NOT SAVE FILES !!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!!!!!! YOU'VE BEEN WARNNED !!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		}
		//clog<<"writing table:"<<myname<<endl;
		//lthresh=l_min_vec[bin_index]/1.5;
		g.set(bin_index,bessel_order);

	}
	fhandler.close();
	// make table of myself and save or load existing one
	//clog<<"LLOW="<<LLOW<<" LHIGH="<<LHIGH<<" NLBINS="<<NLBINS<<endl;
	clog<<myname<<endl;
	loadTable(LLOW,LHIGH,NLBINS,true);
	extrapolationOff();
}


number BandPower_W::get(number ell1)
{
	ell=ell1;
	if(noApodise)
	{
		//clog<<"non apodised values:"<<endl;
		return W_noApodise(ell);
	}
	else
	{
		//clog<<"apodised values"<<endl;
		number result=0.;
		//set integration from thetamin and thetamin+Delta_theta to be different
		//number tmin=thetamin+exp(Delta_x);
		//result=gaussianIntegrate_gsl(*this,thetamin,tmin,100);
		//clog<<"tmin="<<tmin/arcmin<<" result="<<result<<endl;
		//exit(1);
		determine_integration_limits();
		for(unsigned int i=0; (i+1)<integ_limits.size(); i++)
		{
			number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
			result+=res;
		}
		//number tmax=thetamax-exp(Delta_x);
		//result+=gaussianIntegrate_gsl(*this,tmax,thetamax,100);
		return result;
	}
}

int BandPower_W::show_bin_index()
{
	return bin_index;
}

void BandPower_W::determine_integration_limits()
{
	const int Nbins = 10000;
	integ_limits.clear();
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
		table.load(1,i,gsl_sf_bessel_Jn(bessel_order,ell*table.get(0,i))*g.value(table.get(0,i)));
	}
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
}


void BandPower_W::determine_integration_limits_W_noAp()
{
	const int Nbins = 1000;
	integ_limits.clear();
	vector<number> integ_y;
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(l_min_vec[bin_index])+log(l_max_vec[bin_index]/l_min_vec[bin_index])/(Nbins-1.)*i));
		table.load(1,i,G_mu(ell,table.get(0,i),bessel_order));
	}
	integ_limits.push_back(l_min_vec[bin_index]);
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
	integ_limits.push_back(l_max_vec[bin_index]);

	// integ_y.push_back(G_mu(ell,l_max_vec[bin_index],bessel_order));
	// matrix integ_limits_mat(2,integ_limits.size());
	// for(int i=0;i<integ_limits.size(); i++)
	// {
	// 	integ_limits_mat.load(0,i,integ_limits[i]);
	// 	integ_limits_mat.load(1,i,integ_y[i]);
	// }
	//integ_limits_mat.printOut((string("integ_limits_noAp")+toString(bessel_order)+string("_")+toString(ell,2)+string(".ascii")).c_str(),5);
	//table.printOut((string("table_G")+toString(bessel_order)+string("_")+toString(ell,2)+string(".ascii")).c_str(),5);
}


//If no apodisation is set, use this for the integration. 
number BandPower_W::W_noApodise(number ell1)
{
	ell=ell1;
	noApodise=true;
	// if (ell<lthresh)
	// 	return gaussianIntegrate_gsl(*this,LLOW,lthresh,100);

	determine_integration_limits_W_noAp();
	number result=0.;
	for(unsigned int i=0; (i+1)<integ_limits.size(); i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
		result+=res;
	}
	return result;
}

//T(theta) apodises the band powers to reduce oscillations around the edges
//Here x=ln(\theta) and Delta_x is the log width of the apodisation.
//x_l and x_u never go beyond the thetamin and thetamax range given by the user.
number BandPower_W::Apodise(number theta)
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
