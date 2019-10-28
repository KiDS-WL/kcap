#include "COSEBIs.h"

COSEBIs::COSEBIs(){}

COSEBIs::~COSEBIs(){}

COSEBIs::COSEBIs(int nMaximum,number thetamin, number thetamax, int nPairs
	,string WnFolderName,string TnFolderName,string OutputTnFolderName)
{
	initialize(nMaximum,thetamin, thetamax, nPairs, WnFolderName, TnFolderName,OutputTnFolderName);
}

void COSEBIs::initialize(int nMaximum,number thetamin, number thetamax, int nPairs
	,string WnFolderName1,string TnFolderName1,string OutputTnFolderName1)
{
	clog<<"in COSEBIs initialize"<<endl;
	DEn_calculated=false;
	Cov_on=false;
	BCov_on=false;
	WnSet=false;
	TnSet=false;
    EnInteglimitSet=false;
	TpmNotDone=true;
	WnFolderName=WnFolderName1;
	TnFolderName=TnFolderName1;
	OutputTnFolderName=OutputTnFolderName1;

	setEparam(nMaximum,thetamin*arcmin,thetamax*arcmin);
	setZbins(nPairs);

}

void COSEBIs::setZbins(int nPairs1)
{
	nPairs=nPairs1;
	nBins=int((-1+sqrt(1+8*nPairs))/2);
}

void COSEBIs::setEparam(int nMaximum1,number thetamin1, number thetamax1)
{
	clog<<"setting En parameters in COSEBIs"<<endl;
	DEn_calculated=false;
	Cov_on=false;
	WnSet=false;
	TnSet=false;
	TpmNotDone=true;
    EnInteglimitSet=false;
	nMaximum=nMaximum1;
	clog<<"nMaximum="<<nMaximum<<endl;
	thetamin=thetamin1;
	clog<<"theta_min="<<thetamin<<endl;
	thetamax=thetamax1;
	clog<<"theta_max="<<thetamax<<endl;
	LHIGH=high*20./thetamax;
}

void COSEBIs::setWns(int nMaximum)
{
	if(!WnSet)
	{
        EnInteglimitSet=false;
		clog<<"Wn not set setting now:"<<endl;
		WnLog WnL(thetamin,thetamax,nMaximum,TnFolderName,WnFolderName);
		Wn_vec.clear();
		for(int n=0; n<nMaximum; n++)
			Wn_vec.push_back(WnL);
		for(int n=0; n<nMaximum; n++)
			Wn_vec[n].set(n+1);
		WnSet=true;
	}
}

void COSEBIs::setTs(int nMaximum)
{
	if(!TnSet)
	{
		clog<<"setting Tpm in COSEBIs"<<endl;
		TpnLog Tpn(thetamin,thetamax,nMaximum,TnFolderName,OutputTnFolderName);
		for(int n=0; n<nMaximum; n++)
			Tpn_vec.push_back(Tpn);
		for(int n=0; n<nMaximum; n++)
		{
			Tpn_vec[n].PrepareLookupTables(n+1);

		}
		TnSet=true;
	}
}



void COSEBIs::setNoise(number A1,vector<number> sigma_e_vec1,vector<number> nBar_vec)
{
	clog<<"setting noise in COSEBIs"<<endl;
	A=A1;
	noise_vec.clear();
	sigma_e_vec=sigma_e_vec1;
	for(int bin=0; bin<sigma_e_vec.size(); bin++)
	{
		noise_vec.push_back(sigma_e_vec[bin]*sigma_e_vec[bin]/(2.*nBar_vec[bin]));
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
	}
}


void COSEBIs::setNoiseToZero()
{
	clog<<"setting noise to zero in COSEBIs, noise_vec.size()="<<noise_vec.size()<<endl;
	for(int bin=0; bin<noise_vec.size(); bin++)
	{
		noise_vec[bin]=0.;
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
	}
}

// void COSEBIs::setNoise_vec(vector<number> noise_vec1)
// {
// 	noise_vec=noise_vec1;
// 	clog<<"noise_vec is read: "<<noise_vec.size()<<endl;
// 	// exit(1);
// 	for(int bin=0; bin<noise_vec.size(); bin++)
// 		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
// 	clog<<"test"<<endl;
// }

void COSEBIs::setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM)
{
	clog<<"setting Ksi in COSEBIs"<<endl;
	Ksi_m_vec.clear();
	Ksi_p_vec.clear();
	if(InputKsiM.size()== nPairs && InputKsiP.size()== nPairs)
	{
		clog<<"Don't panic, number of redshift bin pairs matches the InputKsi vectors"<<endl;
	}
	else
	{
		clog<<"Panic! number of redshift bin pairs DOES NOT match the InputKsi vectors, exiting now ..."<<endl;
		exit(1);
	}
	for(int r=0;r<nPairs;r++)
	{
		Ksi_p_vec.push_back(function_cosebis());
		Ksi_m_vec.push_back(function_cosebis());
	}
	for(int r=0;r<nPairs;r++)
	{
		Ksi_p_vec[r].loadWithValues(theta,InputKsiP[r],true);
		Ksi_m_vec[r].loadWithValues(theta,InputKsiM[r],true);
		Ksi_p_vec[r].extrapolationOff();
		Ksi_m_vec[r].extrapolationOff();
	}
}

void COSEBIs::setPower(vector<number> log_ell,vector<vector<number> > InputPower)
{
	if(InputPower.size()== nPairs)
	{
		//clog<<"Don't panic, number of redshift bin pairs matches the Input Power vector"<<endl;
	}
	else
	{
		clog<<"Panic! number of redshift bin pairs DOES NOT match the Input Power vector, exiting now ..."<<endl;
		exit(1);
	}
	
	powerspectrum_vec.clear();
	
	for(int r=0; r<InputPower.size(); r++)
		powerspectrum_vec.push_back(function_cosebis());
	for(int r=0; r<InputPower.size(); r++)
	{
		powerspectrum_vec[r].loadWithValues(log_ell,InputPower[r],true);
		powerspectrum_vec[r].extrapolationOff();
	}
}


number COSEBIs::integrant(number l)
{
	if(BCov_on) //integrant for Covariance of Bn
	{
		number integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l);
		return integ;
	}
	else if(Cov_on) //integrant for En covarince
	{
		counter++;
		number integ;
		number power1=powerspectrum_vec[rp1].value(l)+delta1noise;
		number power2=powerspectrum_vec[rp2].value(l)+delta2noise;
		number power3=powerspectrum_vec[rp3].value(l)+delta3noise;
		number power4=powerspectrum_vec[rp4].value(l)+delta4noise;
		number powers=power1*power2+power3*power4;
		integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
		return integ;
	}
	else if(Real) // real space integrant for COSEBIs
	{
		if (realPlus) //T_plus integrant
		{
			number theta=l;
			number Tp=Tpn_vec[nT].TpValue(theta);
			number Ksip=Ksi_p_vec[redshiftPair].value(theta);
			number integ=theta*Tp*Ksip;
			return integ;
		}
		else // T_minus integrant
		{
			number theta=l;
			number Tm=Tpn_vec[nT].TnValue(theta);
 			number integ=theta*Tm*Ksi_m_vec[redshiftPair].value(theta);
			return integ;
		}
	}
	else // Fourier space integrant for COSEBIs
	{
		number power=powerspectrum_vec[redshiftPair].value(l);
		number Wn=Wn_vec[nW].value(l);
		number integ=power*Wn*l;
		return integ;
	}
}


number COSEBIs::ReturnPower(number ell,int rPair)
{
	return powerspectrum_vec[rPair].value(ell);
}


matrix COSEBIs::returnIntegrandForPowerCase(vector <number> ell_vec, int n, int redshift)
{
	matrix integ_mat(2, ell_vec.size());
	BCov_on=false;
	Cov_on=false;
	Real=false;
	realPlus=false;
	setWns(nMaximum);
	nW=n-1;
	redshiftPair=redshift;
	for(int i=0; i<ell_vec.size();i++)
	{
		integ_mat.load(0,i,ell_vec[i]);
		integ_mat.load(1,i,integrant(ell_vec[i]));
	}
	return integ_mat;
}

number COSEBIs::delta(int i1, int j1)
{
	return i1==j1? 1.: 0.;
}

void COSEBIs::determine_integration_limits_En()
{
    if(!EnInteglimitSet)
    {
        const int Nbins = 1000000;
        // free old list
        integ_limits_vec.clear();
        for(int n=0; n<nMaximum; n++)
        {
            integ_limits.clear();
            // make table of integrant values (Wn's only) on a very fine grid
            matrix table(2,Nbins);
            lthresh =  2.*pi/thetamax*(n+1.)/12.;
     		
            for(int i=0;i<Nbins;i++)
            {
                table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
                table.load(1,i,Wn_vec[n].value(table.get(0,i)));
            }
            integ_limits.push_back(lthresh);
            
            for(int i=1;i<Nbins-1;i++)
                if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
                    || (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
                integ_limits.push_back(table.get(0,i));
            integ_limits.push_back(LHIGH);
            integ_limits_vec.push_back(integ_limits);
        }
        EnInteglimitSet=true;
    }
}


number COSEBIs::valueEn(int n)
{
	Cov_on=false;
	Real=false;
	lthresh =  2.*pi/thetamax*n/12.;
	number result= 0.;
	nW=n-1;
	if(LLOW<lthresh)
		result= gaussianIntegrate_gsl(*this,LLOW,lthresh,100);	

	for(unsigned int i=0; (i+1)<integ_limits_vec[nW].size(); i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits_vec[nW][i],integ_limits_vec[nW][i+1],20);
		result+=res;
	}

	return result/2./pi;
}


matrix COSEBIs::calEn()
{
	matrix En(nMaximum*nPairs);
	setWns(nMaximum);
    if(!EnInteglimitSet)
        determine_integration_limits_En();

	for(int r=0; r<nPairs; r++)
	{
		redshiftPair=r;
		for(int n1=nMaximum*r,m=1 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
				En.load(n1,valueEn(m));
		}
	}
	return En;
}

void COSEBIs::determine_integration_limits_En_Tp()
{
	const int Nbins = 1000000;
	// free old list
	integ_limitsTp.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(thetamin)+log(thetamax/thetamin)/(Nbins-1.)*i));
		table.load(1,i,Tpn_vec[nT].TpValue(table.get(0,i)));
	}
	integ_limitsTp.push_back(thetamin);
	for(int i=1;i<Nbins-1;i++)
	{
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		{
			integ_limitsTp.push_back(table.get(0,i));
		}
	}
	integ_limitsTp.push_back(thetamax);
}

/// min max finder
void COSEBIs::determine_integration_limits_En_Tm()
{
	const int Nbins = 100000;
	// free old list
	integ_limitsTm.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(thetamin)+log(thetamax/thetamin)/(Nbins-1.)*i));
		table.load(1,i,Tpn_vec[nT].TnValue(table.get(0,i)));
	}
	integ_limitsTm.push_back(thetamin);
	for(int i=1;i<Nbins-1;i++)//2->1
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		integ_limitsTm.push_back(table.get(0,i));
	integ_limitsTm.push_back(thetamax);
}


matrix COSEBIs::valueEn2PCFs(int n)
{
	Cov_on=false;
	Real=true;
	number IntPlus= 0.;
	number IntMinus= 0.;
	nT=n-1;
	determine_integration_limits_En_Tp();
	determine_integration_limits_En_Tm();
	for(unsigned int i=0;(i+1)<integ_limitsTp.size();i++)
	{
		realPlus=true;
		number IntP=gaussianIntegrate_gsl(*this,integ_limitsTp[i],integ_limitsTp[i+1],10);	
		IntPlus+=IntP;		
	}
	for(unsigned int i=0;(i+1)<integ_limitsTm.size();i++)
	{
		realPlus=false;
		number IntM=gaussianIntegrate_gsl(*this,integ_limitsTm[i],integ_limitsTm[i+1],10);
		IntMinus+=IntM;
	}
	matrix E3(3);
	E3.load(0,(IntPlus+IntMinus)/2.);
	E3.load(1,IntPlus);
	E3.load(2,IntMinus);
	return E3;
}


matrix COSEBIs::calEn2PCFs()
{
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int bin1=0;bin1<nBins;bin1++)
	{
		for(int bin2=bin1;bin2<nBins;bin2++)
		{
			int m=0;
			int p1=calP(nBins,bin1,bin2);
			redshiftPair=p1;
			for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
			{
				clog<<"loading E n="<<n1+1<<endl;
				matrix E3=valueEn2PCFs(m);
				En.load(0,n1,m);
				///E3: En, the plus integral, the minus integral
				En.load(1,n1,E3.get(0));
				En.load(2,n1,E3.get(1));
				En.load(3,n1,E3.get(2));
 				En.load(4,n1,0.5*(E3.get(1)-E3.get(2)));
			}
		}
	}
	return En;
}

matrix COSEBIs::valueEn2PCFsKsiInput(matrix& Ksi_mat,int m,number h)
{
	number IntPlus= 0.;
	number IntMinus= 0.;
	int nTheta=Ksi_mat.rows;
	number intp,intm;
	intp=(Ksi_mat.get(0,0)*Ksi_mat.get(1,0)*Tpm_mat_vec[m].get(0,0))/2.;
	IntPlus=intp;
	intp=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(1,nTheta-1)*Tpm_mat_vec[m].get(0,nTheta-1))/2.;
	IntPlus+=intp;
	intm=(Ksi_mat.get(0,0)*Ksi_mat.get(2,0)*Tpm_mat_vec[m].get(1,0))/2.;
	IntMinus=intm;
	IntMinus+=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(2,nTheta-1)*Tpm_mat_vec[m].get(1,nTheta-1))/2.;
	for(int itheta=1; itheta<(nTheta-1); itheta++)
	{
		intp=Ksi_mat.get(0,itheta)*Ksi_mat.get(1,itheta)*Tpm_mat_vec[m].get(0,itheta);
		IntPlus +=intp;
		intm=Ksi_mat.get(0,itheta)*Ksi_mat.get(2,itheta)*Tpm_mat_vec[m].get(1,itheta);
		IntMinus+=intm;
	}
	
	IntPlus*=h;
	IntMinus*=h;
	number IntTotal=(IntPlus+IntMinus)/2.;
	
	matrix E3(3);
	E3.load(0,IntTotal);
 	E3.load(1,IntPlus);
 	E3.load(2,IntMinus);
	return E3;
}

vector<vector<number> > COSEBIs::readKsi(string FileName,int nColumns)
{
	int nRow=0;
	string str;
	number temp;
	vector<vector<number> > Ksi_vecvec;
	ifstream KsiFile((FileName).c_str());
  	clog<<"reading, "<<FileName<<endl;
	if(KsiFile.fail())
	{
		clog<<"error occured during opening: "<<FileName<<endl;
		exit(1);
	}

	string line;
	while(getline(KsiFile,line))
	{
		stringstream stream(line);
		str=stream.peek();
		if (str=="#")
		{
			//nColumns++;
// 			string line;
// 			getline(KsiFile,line);
// 			clog<<line<<endl;
		}
		else
		{
			int cols_temp=0;
      		number temp;
      		vector<number> Ksi_vec;
      		while(stream>>temp)
      		{
        		Ksi_vec.push_back(temp);
        		cols_temp++;
      		}
      		if(cols_temp>0)
      		{
				Ksi_vecvec.push_back(Ksi_vec);
				nRow++;
			}
		}
	}

	KsiFile.close();
	return Ksi_vecvec;
}

vector<int> COSEBIs::FindMinMaxTheta(vector<vector<number> > Ksi_vec)
{
	int row=0;
	vector<int> index(2);
	clog.precision(10);
	number h=Ksi_vec[1][0]-Ksi_vec[0][0];

	number tmin=thetamin/arcmin;
	number tmax=thetamax/arcmin;
	if(Ksi_vec[0][0]-h>tmin || Ksi_vec[Ksi_vec.size()-1][0]+h<tmax)
	{
		clog<<"the theta range is not sufficient:"
			<<Ksi_vec[0][0]-h<<"   "<<Ksi_vec[Ksi_vec.size()-1][0]+h<<endl;
		exit(1);
	}

	int ntmin=ceil((tmin-Ksi_vec[0][0])/h);
	int ntmax=Ksi_vec.size()-1;

	while(tmin>Ksi_vec[ntmin][0])
		ntmin++;

	while(tmax<Ksi_vec[ntmax][0])
		ntmax--;

	index[0]=ntmin;
	index[1]=ntmax;
	return index;
}

void COSEBIs::FindTnKsiInput(matrix theta_mat)
{

	if(TpmNotDone)
	{

		Tpm_mat_vec.clear();
		int nTheta=theta_mat.size();

		for(int n=0; n<nMaximum; n++)
		{
			matrix Tpm_mat(2,nTheta);
			for(int t=0; t<nTheta; t++)
			{
				number theta=theta_mat.get(t);
				number Tp=Tpn_vec[n].TpValue(theta);
				number Tm=Tpn_vec[n].TnValue(theta);

				Tpm_mat.load(0,t,Tp);
				Tpm_mat.load(1,t,Tm);
			}

			Tpm_mat_vec.push_back(Tpm_mat);
		}
	}
	TpmNotDone=false;
}


matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,int nColumns)
{
 	clog<<"in calE2PCFsFromInputKsi in COSEBIs"<<endl;

 	//number of redshift pairs
	nPairs=FileName.size();
	clog<<"nPairs="<<nPairs<<endl;
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]);
		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}

matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,vector<string> corrFile,int nColumns)
{
	clog<<"in calE2PCFsFromInputKsi in COSEBIs with corrections"<<endl;
	for(unsigned int i=0; i<FileName.size(); i++)
	{
		clog<<"FileName is:"<<FileName[i]<<endl;
		clog<<"Correction FileName is:"<<corrFile[i]<<endl;
	}


	nPairs=FileName.size();
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
		vector<vector<number> > Corr_vec=readKsi(corrFile[r],nColumns);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
 			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[i][1]);//xi_+
 			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[i][1]);//xi_-
		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);

		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}

matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,vector<number> Corr_vec,int nColumns)
{
	clog<<"in calE2PCFsFromInputKsi in COSEBIs with correction_vector"<<endl;
	for(unsigned int i=0; i<FileName.size(); i++)
	{
		clog<<"FileName is:"<<FileName[i]<<endl;
	}
	nPairs=FileName.size();
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[r]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[r]);

		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}


//NOTE: the value of LHIGH is important for the off-diagonals. 
// For better precision use a bigger high
void COSEBIs::determine_integration_limitsCov()
{
/* Idea: find a possibly complete list of consecutive local 
minima/maxima of oscillating integrant and integrate between them
*/
	const int Nbins = 1000000;
	// free old list
	integ_limits.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	number lthresh=pi/thetamax/2.;
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
		table.load(1,i,Wn_vec[nW].value(table.get(0,i))*Wn_vec[mW].value(table.get(0,i)));
	}
// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
	integ_limits.push_back(lthresh);
	for(int i=1;i<Nbins-1;i++)//2->1
		if ((table.get(1,i-1)<table.get(1,i)&& table.get(1,i+1)<table.get(1,i))
		|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		      integ_limits.push_back(table.get(0,i));
	integ_limits.push_back(LHIGH);
}

number COSEBIs::valueCov(int n1,int m1)
{
	Cov_on=true;
	nW=n1-1;
	mW=m1-1;
	number lthresh=pi/thetamax/2.;
	determine_integration_limitsCov();
	number result= gaussianIntegrate_gsl(*this,LLOW,lthresh,20);	
	for(unsigned int i=0;(i+1)<integ_limits.size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
		result+=res;
	}
	Cov_on=false;
	return result/A/2./pi;
}

//lthresh integration doesn't need to be separated, this was done as a test. The results change very little.
number COSEBIs::valueBCov(int n1,int m1)
{
	BCov_on=true;
	nW=n1-1;
	mW=m1-1;
	number lthresh=pi/thetamax/2.;
	determine_integration_limitsCov();
	number result= gaussianIntegrate_gsl(*this,LLOW,lthresh,20);	
	for(unsigned int i=0;(i+1)<integ_limits.size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
		result+=res;
	}
	BCov_on=false;
	return result/A/2./pi;
}


void COSEBIs::readNpairs(vector<string> FileName_vec,int nColumns)
{
	clog<<"in readNpairs in COSEBIs"<<endl;
	clog<<"FileName_vec.size is:"<<FileName_vec.size()<<endl;
	nPairs=FileName_vec.size();
	Npair_mat_vec.clear();
	for(int r=0; r<nPairs; r++)
	{
		vector<vector<number> > Ksi_vec=readKsi(FileName_vec[r],nColumns);
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix Npair_mat(2,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			Npair_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Npair_mat.load(1,i-index[0],Ksi_vec[i][7]);
		}
		Npair_mat_vec.push_back(Npair_mat);
		clog<<"Npair_mat[0]="<<Npair_mat.get(0,0)<<"   Npair_mat[1]="<<Npair_mat.get(1,0)<<endl;
	}
}


number COSEBIs::valueNoiseCov_fromInputNpair(int n1,int m1, int np,int bin1,int bin2)
{
	number result=0.;
	int nT=n1-1;
	int mT=m1-1;
	number deltaTheta=(Npair_mat_vec[np].get(0,2)-Npair_mat_vec[np].get(0,1));

	for(int i=0; i<Npair_mat_vec[np].rows; i++)
	{
		number theta=Npair_mat_vec[np].get(0,i);
		number Npairs=Npair_mat_vec[np].get(1,i);
		number Tp_m=Tpn_vec[mT].TpValue(theta);
		number Tm_m=Tpn_vec[mT].TnValue(theta);
		number Tp_n=Tpn_vec[nT].TpValue(theta);
		number Tm_n=Tpn_vec[nT].TnValue(theta);
		number Npairs_th=2.*deltaTheta*theta*pi*A*neff_vec[bin1]*neff_vec[bin2];
		
		number res=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs;
		
// 		number res_B=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npair_B;
		// 		number res_th=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs_th;
		//    		cout<<theta<<"   "<<Npairs<<"    "<<Npairs_th<<endl;
		// 		if( (i % 1000)==0. )
		// 		{
		// 			clog<<"theta="<<theta<<"  Npairs="<<Npairs<<",  Npairs_th="<<Npairs_th<<endl;
		// 			clog<<"res="<<res<<endl;
		// 		}
		//  		result+=res;
		result+=res;
		//  		result+=res_th;
	}

	return result*deltaTheta*deltaTheta;
}

matrix COSEBIs::calNoiseCov_fromInputNpair()
{
	clog<<"calculating the noise covariance in COSEBIs from input npair"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);

	clog<<"Noise Cov file failed"<<endl;
	setTs(nMaximum);
	for(int n=1; n<nMaximum+1; n++)
	{
		for(int m=n; m<nMaximum+1;m++)
		{
			clog.precision(10);

			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							// 							delta1noise=delta(bin1,bin3)*SigmaE_vec[bin1]*SigmaE_vec[bin1];
							// 							delta2noise=delta(bin2,bin4)*SigmaE_vec[bin2]*SigmaE_vec[bin2];
							// 							delta3noise=delta(bin1,bin4)*SigmaE_vec[bin1]*SigmaE_vec[bin1];
							// 							delta4noise=delta(bin2,bin3)*SigmaE_vec[bin2]*SigmaE_vec[bin2];
							// 							cout<<"delta1noise="<<delta1noise<<endl;
							// 							cout<<"delta2noise="<<delta2noise<<endl;
							// 							cout<<"delta3noise="<<delta3noise<<endl;
							// 							cout<<"delta4noise="<<delta4noise<<endl;
							//the pair considered for the Eparam_vec
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							int i=nMaximum*p1+n-1;
							int j=nMaximum*p2+m-1;
							number deltas=delta(bin1,bin3)*delta(bin2,bin4)+delta(bin1,bin4)*delta(bin2,bin3);
							///Why this? This is because of the way Npair is accounted for in Athena. They are counted twice when i=j
							if(bin1==bin2)
								deltas=deltas*0.5;
							
							number CovNM=0.;
							if (deltas>0)
								CovNM=valueNoiseCov_fromInputNpair(n,m,p1,bin1,bin2);

							number CovNoise=CovNM*deltas*
								sigma_e_vec[bin1]*sigma_e_vec[bin1]*sigma_e_vec[bin1]*sigma_e_vec[bin1]/8.;
							CMT.load(i,j,CovNoise);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}
	return CMT;
}


bool COSEBIs::checkPower()
{
	if(powerspectrum_vec.size())
		return true;
	return false;
}

matrix COSEBIs::calCov()
{
	clog<<"calculating the covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	//clog<<"calculating the covariance in COSEBIs"<<endl;
	if(checkPower())
	{
		clog<<"power spectra are set, nBins="<<nBins<<endl;
		for(int bin1=0; bin1<nBins; bin1++)
		{
			for(int bin2=bin1; bin2<nBins; bin2++)
			{
				for(int bin3=bin1; bin3<nBins; bin3++)
				{
					for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
					{
						//clog<<"noise_vec="<<noise_vec[bin1]<<"  "<<noise_vec[bin2]<<endl;
						int n,m=0;
						rp1=calP(nBins,bin1,bin3);
						rp2=calP(nBins,bin2,bin4);
						rp3=calP(nBins,bin1,bin4);
						rp4=calP(nBins,bin2,bin3);
						delta1noise=delta(bin1,bin3)*noise_vec[bin1];
						delta2noise=delta(bin2,bin4)*noise_vec[bin2];
						delta3noise=delta(bin1,bin4)*noise_vec[bin1];
						delta4noise=delta(bin2,bin3)*noise_vec[bin2];
						//the pair considered for the Eparam_vec
						int p1=calP(nBins,bin1,bin2);
						int p2=calP(nBins,bin3,bin4);
						for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
						{
							for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
							{
								CMT.load(i,j,valueCov(n,m));
								CMT.load(j,i,CMT.get(i,j));
								CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
								CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
							}
						}
					}
				}
			}
		}
	}
	else
	{
		clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		clog<<"!!!!!WARNING power spectrum is not set. The results are not reliable!!!!!!!!!!!!!!"<<endl;
		clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	}
	return CMT;
}


matrix COSEBIs::calCovForSigma_m(number sigma_m)
{
	clog<<"calculating the covariance for sigma_m in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	matrix En_mat=calEn();
	for(int i=0;i<nMaximum*(nPairs);i++)
	{
		for(int j=0; j<nMaximum*(nPairs); j++)
		{
			CMT.load(i,j,En_mat.get(i)*En_mat.get(j));
		}
	}

	return CMT*4.*sigma_m*sigma_m;
}



matrix COSEBIs::calBCov()
{
	clog<<"calculating the B-mode covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	clog<<" BCov file failed"<<endl;
	for(int n=1; n<nMaximum+1; n++)
	{
		for(int m=n; m<nMaximum+1;m++)
		{
			//clog<<"n="<<n<<"  m="<<m<<endl;
			number CovNM=valueBCov(n,m);
			clog.precision(10);
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];

							//the pair considered for the Eparam_vec
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							int i=nMaximum*p1+n-1;
							int j=nMaximum*p2+m-1;
// 								cout<<"i="<<i<<"  j="<<j<<endl;
// 								bool writeInteg=false;
							number CovB=CovNM*(delta1noise*delta2noise+delta3noise*delta4noise);
							CMT.load(i,j,CovB);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}
	return CMT;
}
///////////////////////////covariance functions ended//////////////////////////////////////////


int COSEBIs::calP(int nBins,int fbin,int sbin)
{
	if(fbin>sbin)
	{
		int swap=fbin;
		fbin=sbin;
		sbin=swap;
	}
	int p=fbin*nBins;

	for(int i=0; i<fbin; i++)
	  	p-=i;
	return p+sbin-fbin;
}


