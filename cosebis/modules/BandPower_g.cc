#include "BandPower_g.h"

BandPower_g::BandPower_g()
{
	thetamin=1.;
	thetamax=1.;
	LLOW=1.;
	LHIGH=1e6;
	FolderName="./cosebis/BandPower/";
	FileName="g";
}

//constructor 
BandPower_g::BandPower_g(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,number LLOW, number LHIGH,
		string FolderName,string FileName)
{
  clog<<"in BandPower_g.........."<<endl;
  initialize(thetamin,thetamax,Response_type,l_min_vec,l_max_vec,LLOW,LHIGH,FolderName,FileName);
  setBandPower_gName(FolderName,FileName);
  setTheta(thetamin,thetamax);
  Set_Response_function(Response_type,l_min_vec,l_max_vec);
  Analytic=false;
}

void BandPower_g::initialize(number thetamin,number thetamax, string Response_type,
		vector<number> l_min_vec,vector<number> l_max_vec,number LLOW, number LHIGH,
		string FolderName,string FileName)
{
  clog<<"in BandPower_g initialize"<<endl;
  setTheta(thetamin,thetamax);
  setBandPower_gName(FolderName,FileName);
  Set_Response_function(Response_type,l_min_vec,l_max_vec);
  Analytic=false;
  set_l_min_max(LLOW,LHIGH);
}

BandPower_g::~BandPower_g(){}

//set thetamin and thetamax
void BandPower_g::setTheta(number thetamin1,number thetamax1)
{
	//clog<<"setting thetamin,thetamax"<<endl;
	thetamin=thetamin1;
	thetamax=thetamax1;
}

//set bandpower folder and file name
void BandPower_g::setBandPower_gName(string FolderName1,string FileName1)
{
	FolderName=FolderName1;
	FileName=FileName1;
}

//set the parameters for the response function S^i(\ell)
void BandPower_g::Set_Response_function(string type, vector<number> l_min_vec1,
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
}

//set the minimum and maximum ell used overall
void BandPower_g::set_l_min_max(number LLOW1,number LHIGH1)
{
	LLOW=LLOW1;
	LHIGH=LHIGH1;
}

//set analytic g_\mu, currently only know the analytic solution for the tophat response function
void BandPower_g::setAnalytic(bool Analytic1)
{
	Analytic=Analytic1;
	if(Analytic)
	{
		//clog<<"setting analytic solution for g_pm"<<endl;
		if(Response_function_type=="tophat")
		{
			//clog<<"For the top hat response function S_i(l)"<<endl;
		}
		else
		{
			clog<<"Don't know the analytic solution for "<<Response_function_type<<" function, setting to false"<<endl;
			Analytic=false;
		}
	}
}

//the response function value, S^i(\ell)
//currently only top hat is defined
number BandPower_g::S_response(int bin_index,number l)
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

//write the reponse function for band power bin_index=i 
//S^i(\ell) for ell between the total maximum and minum ell value
void BandPower_g::writeResponse(int bin_index)
{
	matrix S_mat(2,10000);
	for(int i=0; i<10000; i++)
	{
		number ell=exp(log(LLOW)+log(LHIGH/LLOW)/(10000-1)*i);
		S_mat.load(0,i,ell);
		S_mat.load(1,i,S_response(bin_index,ell));
	}
	bin_index++;
	//clog<<"writting the response function ......"<<endl;
	//clog<<"bin_index="<<bin_index<<endl;
	S_mat.printOut((FolderName+string("/S_")+toString(bin_index)+string("_")
		+toString(nBands)+string("_")+toString(l_min_vec[0],2)
		+string("-")+toString(l_max_vec[nBands-1],2)+string(".ascii")).c_str(),8);
}


//here we set the g function
void BandPower_g::set(int bin_index1,int bessel_order1)
{
	//clog<<"set bin_index and bessel_order for BandPower_g"<<endl;
	bin_index=bin_index1;
	bessel_order=bessel_order1;
	// is there a table on disk?
	string pm_string=toString(bessel_order);

	string myname =FolderName+string("/")+FileName+string("_")+pm_string
			+string("_")+toString(bin_index+1)+string("_")
		 	+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
	ifstream fhandler((myname+string(".table")).c_str());
	// NO?
	if (fhandler.fail())
	{
		if(!CheckFolderExist(FolderName))
		{
			//clog<<"making the folder for BandPower_g:"<<FolderName<<endl;
			mkdir((FolderName).c_str(), 0777);
		}
		clog<<"writing table:"<<myname<<endl;
		determine_integration_limits();
	}
	fhandler.close();
	// make table of myself and save or load existing one
	clog<<"thetamin="<<thetamin<<" thetamax="<<thetamax<<" nTheta="<<nTheta<<endl;
	makeTable(thetamin,thetamax,nTheta,false);
	extrapolationOff();
}


//The integrant is only used if analytic is set to false
number BandPower_g::integrant(number x)
{
	number integ=x*S_response(bin_index,x/theta)*gsl_sf_bessel_Jn(bessel_order,x);
	return integ;
}

number BandPower_g::get(number theta1)
{
	theta=theta1;
	if(Analytic)
	{
		return Theta_g_tophat(theta);
	}
	else
	{
		number Integral=0.;

		if(Response_function_type=="tophat")
		{
			//clog<<"reponse is tophat"<<endl;
			number l_min=l_min_vec[bin_index];
			number l_max=l_max_vec[bin_index];
			int i_start=0;
			for(i_start=0; integ_limits[i_start]<l_min*theta; i_start++);

			Integral=gaussianIntegrate_gsl(*this,l_min*theta,integ_limits[i_start],10);
			//clog<<"Integral="<<Integral<<endl;
			int i=i_start;	
			while(integ_limits[i]<l_max*theta)
			{
				number Int=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],10);
				Integral+=Int;
				//clog<<"i="<<i<<" Intgeral="<<Integral<<endl;
				i++;

			}
			//number Int=gaussianIntegrate_gsl(*this,l_max*theta,integ_limits[i],10);
		}
		else
		{
			for(int i=0; i<integ_limits.size(); i++)
			{
				number Int=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],10);
				Integral+=Int;
			}
		}
		return Integral/theta/theta;
	}
}


int BandPower_g::show_bin_index()
{
	return bin_index;
}

//check the range
void BandPower_g::determine_integration_limits()
{
	const int Nbins = 1000000;
	// free old list
	integ_limits.clear();
	number x_min=LLOW*thetamin;
	number x_max=LHIGH*thetamax;
	if(Response_function_type=="tophat")
	{
		x_min=l_min_vec[0]*thetamin;
		x_max=l_max_vec[nBands-1]*thetamax;
	}
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(x_min)+log(x_max/x_min)/(Nbins-1.)*i));
		table.load(1,i,gsl_sf_bessel_Jn(bessel_order,table.get(0,i)));
	}
	integ_limits.push_back(x_min);
	for(int i=1;i<Nbins-1;i++)
	{
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		{
			integ_limits.push_back(table.get(0,i));
		}
	}
	integ_limits.push_back(x_max);
}

//if analytic then uses this one
number BandPower_g::Theta_g_tophat(number theta)
{
	if(bessel_order==0)
	{
		number ell_theta_u=l_max_vec[bin_index]*theta;
		number ell_theta_l=l_min_vec[bin_index]*theta;
		number result= 1./theta*(
			 l_max_vec[bin_index]*gsl_sf_bessel_Jn(1,ell_theta_u)
			-l_min_vec[bin_index]*gsl_sf_bessel_Jn(1,ell_theta_l));
		return result;
	}
	else if(bessel_order==2)
	{
		number ell_theta_u=l_max_vec[bin_index]*theta;
		number ell_theta_l=l_min_vec[bin_index]*theta;
		number result=-1./theta/theta*(
			 ell_theta_u*gsl_sf_bessel_Jn(1,ell_theta_u)
			-ell_theta_l*gsl_sf_bessel_Jn(1,ell_theta_l)
			+2.*gsl_sf_bessel_Jn(0,ell_theta_u)
			-2.*gsl_sf_bessel_Jn(0,ell_theta_l));
		return result;
	}
	else if(bessel_order==4)
	{
		number ell_theta_u=l_max_vec[bin_index]*theta;
		number ell_theta_l=l_min_vec[bin_index]*theta;
		number G_m_ell_theta_u=((ell_theta_u-8./ell_theta_u)*gsl_sf_bessel_Jn(1,ell_theta_u))-8.*gsl_sf_bessel_Jn(2,ell_theta_u);
		number G_m_ell_theta_l=((ell_theta_l-8./ell_theta_l)*gsl_sf_bessel_Jn(1,ell_theta_l))-8.*gsl_sf_bessel_Jn(2,ell_theta_l);
		number result= 1./theta/theta*(G_m_ell_theta_u-G_m_ell_theta_l);
		return result;
	}
	else
	{
		return 0;
	}
}
