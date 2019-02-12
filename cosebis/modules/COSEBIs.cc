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
	//isFlat=false;
	TpmNotDone=true;
	WnFolderName=WnFolderName1;
	TnFolderName=TnFolderName1;
	OutputTnFolderName=OutputTnFolderName1;
	setEparam(nMaximum,thetamin*arcmin,thetamax*arcmin);
	setZbins(nPairs);
	//setNoise(A,sigmaE,nBar);
}

void COSEBIs::setZbins(int nPairs1)
{
	nPairs=nPairs1;
}

void COSEBIs::setEparam(int nMaximum1,number thetamin1, number thetamax1)
{
	clog<<"setting En parameters in COSEBIs"<<endl;
	DEn_calculated=false;
	Cov_on=false;
	WnSet=false;
	TnSet=false;
	//isFlat=false;
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
		WnLog WnL(thetamin,thetamax,nMaximum,WnFolderName,TnFolderName);
		Wn_vec.clear();
		///these need to be in separate for loops, otherwise a segmentation fault happens!
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



void COSEBIs::setNoise(number A1,number sigmaE1,number nBar1)
{
	clog<<"setting noise in COSEBIs"<<endl;
	A=A1;
	sigmaE=sigmaE1;
	nBar=nBar1;
	number nGal=nBar/nBins;
	noise_vec.clear();
	for(int bin=0; bin<nBins; bin++)
	{
		noise_vec.push_back(sigmaE*sigmaE/2./nGal);
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
	}
}

void COSEBIs::setNoise_vec(vector<number> noise_vec1)
{
	noise_vec=noise_vec1;
	for(int bin=0; bin<nBins; bin++)
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
}

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
	clog<<"in setPower in COSEBIs.cc"<<endl;
	//check if InputPower.size()== nPairs
	if(InputPower.size()== nPairs)
	{
		clog<<"Don't panic, number of redshift bin pairs matches the Input Power vector"<<endl;
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
// 		string powerFileName=string("Power_")+toString(r);
//   		powerspectrum_vec[r].setName(powerFileName.c_str(),function_cosebis::NONAMECOUNTER);
//  		powerspectrum_vec[r].makeTable(0.15,100000.,200,true);
//  		powerspectrum_vec[r].saveTable();
	}
}

matrix COSEBIs::readInputEn(string InputEnfileName)
{
	En_data.restore(InputEnfileName.c_str());
	return En_data;
}


matrix COSEBIs::readInputCovariance(string InputCovfileName)
{
	Cov_mat.restore(InputCovfileName.c_str());
	return Cov_mat;
}

number COSEBIs::CalChiS(matrix En_th,matrix En_data,matrix Cov_mat)
{
// 	clog<<"En_th.rows="<<En_th.rows<<" En_data.rows="<<En_data.rows<<endl;
// 	clog<<"En_th.columns="<<En_th.columns<<" En_data.columns="<<En_data.columns<<endl;
// 	clog<<"Cov_mat.columns="<<Cov_mat.columns<<" Cov_mat.rows="<<Cov_mat.rows<<endl;
	matrix DeltaEn=En_data-En_th;
	matrix iCov=Cov_mat.inverse();
	matrix chiS=DeltaEn.t()*iCov*DeltaEn;
// 	clog<<"chiS.size()="<<chiS.size()<<endl;
// 	clog<<"chiS.get(0)="<<chiS.get(0)<<endl;
	return chiS.get(0);
}



number COSEBIs::integrant(number l)
{
	if(BCov_on)
	{
		number integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l);
		return integ;
	}
	else if(Cov_on)
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
	else if(Real)
	{
		if (realPlus)
		{
			number theta=l;
			number Tp=Tpn_vec[nT].TpValue(theta);
			number Ksip=Ksi_p_vec[redshiftPair].value(theta);
			number integ=theta*Tp*Ksip;
			return integ;
		}
		else
		{
			number theta=l;
			number Tm=Tpn_vec[nT].TnValue(theta);
 			number integ=theta*Tm*Ksi_m_vec[redshiftPair].value(theta);
			return integ;
		}
	}
	else
	{
// 		clog<<"In integrant"<<endl;
		number power=powerspectrum_vec[redshiftPair].value(l);
		number Wn=Wn_vec[nW].value(l);
		number integ=power*Wn*l;
// 		cout<<l<<"    "<<power<<"   "<<Wn<<endl;
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
    // 		clog<<"lthresh="<<lthresh<<" LHIGH="<<LHIGH<<endl;
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
// 	clog<<"In valueEn, lthresh="<<lthresh<<endl;
	nW=n-1;
	number result= gaussianIntegrate_gsl(*this,LLOW,lthresh,100);	
// 	clog<<"result="<<result<<endl;
// 	clog<<"n="<<n<<"  integ_limits.size()="<<integ_limits.size()<<endl;
	for(unsigned int i=0; (i+1)<integ_limits_vec[nW].size(); i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits_vec[nW][i],integ_limits_vec[nW][i+1],20);
		result+=res;
	}
// 	clog<<"result="<<result<<endl;
	return result/2./pi;
}

matrix COSEBIs::calEn()
{
	//clog<<"nPairs="<<nPairs<<endl;
	matrix En(nMaximum*nPairs);
	setWns(nMaximum);
    if(!EnInteglimitSet)
        determine_integration_limits_En();
//     else
//         clog<<"integ limits are already set"<<endl;
	for(int r=0; r<nPairs; r++)
	{
		redshiftPair=r;
		for(int n1=nMaximum*r,m=1 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
				//clog<<"loading E n="<<n1+1<<endl;
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
	// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
	//   clog<<"---------------------------------------------------"<<endl;
	//   clog<<"COSEBIs"<<endl;
	//   clog<<"LLOW="<<LLOW<<endl;
	//   clog<<"LHIGH="<<LHIGH<<endl;
	//   clog<<"---------------------------------------------------"<<endl;
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
	//cout<<"#integlimitsTp_size="<<integ_limitsTp.size()<<endl;
// 	table.printOut((string("integLimitTableTp-")+string("-")+toString(nT+1)+string("_")
// 	+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)
// 	+string(".ascii")).c_str(),10);
	//clog<<"integ_limits En done"<<endl;
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
	// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
	//   clog<<"---------------------------------------------------"<<endl;
	//   clog<<"COSEBIs"<<endl;
	//   clog<<"LLOW="<<LLOW<<endl;
	//   clog<<"LHIGH="<<LHIGH<<endl;
	//   clog<<"---------------------------------------------------"<<endl;
	integ_limitsTm.push_back(thetamin);
	for(int i=1;i<Nbins-1;i++)//2->1
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		integ_limitsTm.push_back(table.get(0,i));
	integ_limitsTm.push_back(thetamax);
	//cout<<"#integlimitsTm_size="<<integ_limitsTm.size()<<endl;
// 	table.printOut((string("integLimitTableTm-")+string("-")+toString(nT+1)+string("_")
// 	+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)
// 	+string(".ascii")).c_str(),10);
	//clog<<"integ_limits En done"<<endl;
}


matrix COSEBIs::valueEn2PCFs(int n)
{
	//cout.precision(5);
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
	clog<<"in calE2PCFs in COSEBIs"<<endl;
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
//  	clog<<"reading, "<<FileName<<endl;
	if(KsiFile.fail())
	{
		clog<<"error occured during opening: "<<FileName<<endl;
		exit(1);
	}
	///lets change this
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
			vector<number> Ksi_vec(nColumns);
			for(int nCol=0; nCol<nColumns; nCol++)
			{
				stream>>temp;
				Ksi_vec[nCol]=temp;
			}
			Ksi_vecvec.push_back(Ksi_vec);
			nRow++;
		}
	}
	KsiFile.close();
	return Ksi_vecvec;
}

vector<int> COSEBIs::FindMinMaxTheta(vector<vector<number> > Ksi_vec)
{
	int row=0;
// 	clog<<"finding min and max theta that are closest to thetamin and thetamax"<<endl;
	vector<int> index(2);
	clog.precision(10);
	number h=Ksi_vec[1][0]-Ksi_vec[0][0];
// 	clog<<"h="<<h<<endl;
	number tmin=thetamin/arcmin;
	number tmax=thetamax/arcmin;
	if(Ksi_vec[0][0]>tmin || Ksi_vec[Ksi_vec.size()-1][0]<tmax)
	{
		clog<<"the theta range is not sufficient:"
			<<Ksi_vec[0][0]<<"   "<<Ksi_vec[Ksi_vec.size()-1][0]<<endl;
		exit(1);
	}
	//clog<<"Ksi_vec[900][0]="<<Ksi_vec[900][0]<<endl;
	//clog<<"Ksi_vec[901][0]="<<Ksi_vec[901][0]<<endl;
	
	//clog<<"tmin="<<tmin<<endl;
	//clog<<"tmax="<<tmax<<endl;
	int ntmin=ceil((tmin-Ksi_vec[0][0])/h);
	int ntmax=Ksi_vec.size()-1;
	//clog<<"at the start ntmin="<<ntmin<<endl;
	//clog<<"tmax-Ksi_vec[0][0]="<<tmax-Ksi_vec[0][0]<<endl;
	//clog<<"at the start ntmax="<<ntmax<<endl;
	while(tmin>Ksi_vec[ntmin][0])
		ntmin++;
	//clog<<"Ksi_vec["<<ntmin<<"]="<<Ksi_vec[ntmin][0]<<endl;
	while(tmax<Ksi_vec[ntmax][0])
		ntmax--;
	//ntmax--;
	//clog<<"ntmax="<<ntmax<<endl;
	//clog<<"Ksi_vec["<<ntmax<<"]="<<Ksi_vec[ntmax][0]<<endl;
	index[0]=ntmin;
	index[1]=ntmax;
	return index;
}

void COSEBIs::FindTnKsiInput(matrix theta_mat)
{
// 	clog<<"in FindTnKsiInput"<<endl;
	if(TpmNotDone)
	{
// 		clog<<"TpmNotDone"<<endl;
		Tpm_mat_vec.clear();
		int nTheta=theta_mat.size();
// 		clog<<"nTheta="<<nTheta<<endl;
		for(int n=0; n<nMaximum; n++)
		{
			matrix Tpm_mat(2,nTheta);
			for(int t=0; t<nTheta; t++)
			{
				number theta=theta_mat.get(t);
				number Tp=Tpn_vec[n].TpValue(theta);
				number Tm=Tpn_vec[n].TnValue(theta);
				//clog<<"theta="<<theta<<"  Tp_n(theta)="<<Tp<<"  Tm_n(theta)="<<Tm<<endl;
				//Tpm_mat.load(0,t,theta);
				Tpm_mat.load(0,t,Tp);
				Tpm_mat.load(1,t,Tm);
			}
			//Tpm_mat.printOut((string("Tpm_mat_")+toString(n+1)+string(".ascii")).c_str(),10);
			Tpm_mat_vec.push_back(Tpm_mat);
		}
	}
	TpmNotDone=false;
}

///check this!
matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,int nColumns)
{
// 	clog<<"in calE2PCFsFromInputKsi in COSEBIs"<<endl;
// 	for(unsigned int i=0; i<FileName.size(); i++)
// 		clog<<"FileName is:"<<FileName[i]<<endl;

	//number deltaTheta=thetamax-thetamin;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	//number of redshift pairs
	nPairs=FileName.size();
// 	clog<<"nPairs="<<nPairs<<endl;
	nBins=(sqrt(8*nPairs+1)-1)/2;
	for(int r=0;r<nPairs;r++)
	{
// 		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
// 		clog<<"index[0]="<<index[0]<<"  index[1]="<<index[1]<<endl;
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]);
		}
		//clog<<"theta1="<<Ksi_mat.get(0,0)/arcmin<<"  theta2="<<Ksi_mat.get(0,1)/arcmin<<endl;
		//clog<<"Ksip(theta1)="<<Ksi_mat.get(1,0)<<"   Ksip(theta2)="<<Ksi_mat.get(1,1)<<endl;
		//clog<<"Ksim(theta1)="<<Ksi_mat.get(2,0)<<"   Ksim(theta2)="<<Ksi_mat.get(2,1)<<endl;
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		//clog<<"h="<<h<<endl;
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

	//number deltaTheta=thetamax-thetamin;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	//number of redshift pairs
	nPairs=FileName.size();
	//clog<<"nPairs="<<nPairs<<endl;
	nBins=(sqrt(8*nPairs+1)-1)/2;
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
		vector<vector<number> > Corr_vec=readKsi(corrFile[r],nColumns);
// 		clog<<Corr_vec[0][0]<<endl;
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		//clog<<"index[0]="<<index[0]<<"  index[1]="<<index[1]<<endl;
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		//matrix Ksi_matnoCorr(3,index[1]-index[0]+1);
		//matrix Corr_mat(2,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
 			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[i][1]);//xi_+
 			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[i][1]);//xi_-

			
			//Ksi_matnoCorr.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			//Ksi_matnoCorr.load(1,i-index[0],Ksi_vec[i][1]);
			//Ksi_matnoCorr.load(2,i-index[0],Ksi_vec[i][2]);
			
			//Corr_mat.load(0,i-index[0],Corr_vec[i][0]*arcmin);
			//Corr_mat.load(1,i-index[0],Corr_vec[i][1]);
		}
// 		Ksi_mat.printOut((folderName+string("/Ksi_mat.ascii")).c_str(),10);
// 		exit(1);
		//Ksi_matnoCorr.printOut((folderName+string("/Ksi_matnoCorr.ascii")).c_str(),10);
		//Corr_mat.printOut((folderName+string("/Corr_mat.ascii")).c_str(),10);
		
		//exit(1);		
		//clog<<"theta1="<<Ksi_mat.get(0,0)/arcmin<<"  theta2="<<Ksi_mat.get(0,1)/arcmin<<endl;
		//clog<<"Ksip(theta1)="<<Ksi_mat.get(1,0)<<"   Ksip(theta2)="<<Ksi_mat.get(1,1)<<endl;
		//clog<<"Ksim(theta1)="<<Ksi_mat.get(2,0)<<"   Ksim(theta2)="<<Ksi_mat.get(2,1)<<endl;
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		//clog<<"h="<<h<<endl;
		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
// 			clog<<E3.get(0)<<'\t'<<E3.get(1)<<'\t'<<E3.get(2)<<'\t'<<(E3.get(1)-E3.get(2))/2.<<endl;
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

	//number deltaTheta=thetamax-thetamin;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	//number of redshift pairs
	nPairs=FileName.size();
	//clog<<"nPairs="<<nPairs<<endl;
	nBins=(sqrt(8*nPairs+1)-1)/2;
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
// 		vector<vector<number> > Corr_vec=readKsi(corrFile[r],nColumns);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		//clog<<"index[0]="<<index[0]<<"  index[1]="<<index[1]<<endl;
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		//matrix Ksi_matnoCorr(3,index[1]-index[0]+1);
		//matrix Corr_mat(2,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[r]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[r]);
			
			//Ksi_matnoCorr.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			//Ksi_matnoCorr.load(1,i-index[0],Ksi_vec[i][1]);
			//Ksi_matnoCorr.load(2,i-index[0],Ksi_vec[i][2]);
			
			//Corr_mat.load(0,i-index[0],Corr_vec[i][0]*arcmin);
			//Corr_mat.load(1,i-index[0],Corr_vec[i][1]);
		}
		//Ksi_mat.printOut((folderName+string("/Ksi_mat.ascii")).c_str(),10);
		//Ksi_matnoCorr.printOut((folderName+string("/Ksi_matnoCorr.ascii")).c_str(),10);
		//Corr_mat.printOut((folderName+string("/Corr_mat.ascii")).c_str(),10);
		
		//exit(1);		
		//clog<<"theta1="<<Ksi_mat.get(0,0)/arcmin<<"  theta2="<<Ksi_mat.get(0,1)/arcmin<<endl;
		//clog<<"Ksip(theta1)="<<Ksi_mat.get(1,0)<<"   Ksip(theta2)="<<Ksi_mat.get(1,1)<<endl;
		//clog<<"Ksim(theta1)="<<Ksi_mat.get(2,0)<<"   Ksim(theta2)="<<Ksi_mat.get(2,1)<<endl;
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		//clog<<"h="<<h<<endl;
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




//NOTE: the value of LHIGH is important for the off-diagonals. For better precision use a bigger high
void COSEBIs::determine_integration_limitsCov()
{
/* Idea: find a possibly complete list of consecutive local minima/maxima of oscillating integrant and integrate between them
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

matrix COSEBIs::calCov()
{
	clog<<"calculating the covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	clog<<" CMT file failed"<<endl;
	//calculates all the powers
	//need to replace
	//calPower();
	for(int bin1=0; bin1<nBins; bin1++)
	{
		for(int bin2=bin1; bin2<nBins; bin2++)
		{
			for(int bin3=bin1; bin3<nBins; bin3++)
			{
				for(int bin4=max(bin2,bin3); bin4<nBins;bin4++)
				{
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
						//cout<<"in COSEBIs calCov, bin1="<<bin1<<"  bin2="<<bin2
						//	<<"  bin3="<<bin3<<"  bin4="<<bin4<<"  p1="<<p1<<"  p2="<<p2<<endl;
						//cout<<"rp1="<<rp1<<" rp2="<<rp2<<" rp3="<<rp3<<" rp4="<<rp4<<endl;
// 						cout<<"noise1="<<delta1noise<<" noise2="<<delta2noise<<" noise3="<<delta3noise
// 							<<"noise4="<<delta4noise<<endl;
					for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
// 								bool writeInteg=false;
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
	return CMT;
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
						for(int bin4=max(bin2,bin3); bin4<nBins;bin4++)
						{
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
// 								cout<<"delta1noise="<<delta1noise<<endl;
// 								cout<<"delta2noise="<<delta2noise<<endl;
// 								cout<<"delta3noise="<<delta3noise<<endl;
// 								cout<<"delta4noise="<<delta4noise<<endl;
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


/*number COSEBIs::trapIntegrant(number theta)
{
	if (realPlus)
	{
		number Tp=Tpn_vec[nT].TpValue(theta);
		number Ksip=Ksi_p_vec[redshiftPair].value(theta);
		if(noisyKsi)
			Ksip+=NoiseKsi/sqrt(theta)*random_matP.get(iRand);
			//number integ=Tp*Tp;
		number integ=theta*Tp*Ksip;
		//if(nT==1)
 		//	cout<<theta<<"\t"<<Ksip<<endl;
		return integ;
	}
	else
	{
		number Tm=Tpn_vec[nT].TnValue(theta);
		number Ksim=Ksi_m_vec[redshiftPair].value(theta);
		if(noisyKsi)
			Ksim+=NoiseKsi/sqrt(theta)*random_matM.get(iRand);
			//number integ=Tm*Tm;
 		number integ=theta*Tm*Ksim;
// 			if(nT==19)
// 				cout<<theta<<"\t"<<integ<<endl;
		return integ;
	}
}*/

/*matrix COSEBIs::valueEn2PCFsTrap(int n,int nBinsTrap,number h)
{
	//cout.precision(5);
	Cov_on=false;
	Real=true;
	number IntPlus= 0.;
	number IntMinus= 0.;
	nT=n-1;
	//cout<<"#n="<<n<<endl;
	realPlus=true;
	iRand=0;
	number integThetaminPlus=trapIntegrant(thetamin);
	iRand=nBinsTrap-1;
	number integThetamaxPlus=trapIntegrant(thetamax);
	//cout<<"#integThetaminPlus="<<integThetaminPlus<<endl;
	//cout<<"#integThetamaxPlus="<<integThetamaxPlus<<endl;
	realPlus=false;
	iRand=0;
	number integThetaminMinus=trapIntegrant(thetamin);
	iRand=nBinsTrap-1;
	number integThetamaxMinus=trapIntegrant(thetamax);
	//cout<<"#integThetaminMinus="<<integThetaminMinus<<endl;
	//cout<<"#integThetamaxMinus="<<integThetamaxMinus<<endl;
	for(int i=1;i<(nBinsTrap-1);i++)
	{
		//clog<<"i="<<i<<endl;
		iRand=i;
		realPlus=true;
		number t=h*number(i)+thetamin;
		number IntP=trapIntegrant(t);
		realPlus=false;
		number IntM=trapIntegrant(t);
		IntMinus+=IntM;
		IntPlus+=IntP;
		//cout<<i<<"\t"<<t<<"\t"<<IntPlus<<"\t"<<IntMinus<<endl;
	}

	IntPlus+=integThetaminPlus/2.;
	IntPlus+=integThetamaxPlus/2.;
	IntMinus+=integThetaminMinus/2.;
	IntMinus+=integThetamaxMinus/2.;
	IntPlus=IntPlus*h;
	IntMinus=IntMinus*h;
	matrix E3(3);
	//E3.load(0,n);
	E3.load(0,(IntPlus+IntMinus)/2.);
	E3.load(1,IntPlus);
	E3.load(2,IntMinus);
	return E3;
}*/


///the noise part only works for no tomography for now and is very simplified
/*vector<matrix> COSEBIs::calEn2PCFsTrap(int nBinsTrapBegin,int nBinsTrapEnd,int stepSize,
								number sigma8,number omega_matter,number omega_lambda,
								number w0,number nindex,number h100,number omega_baryon,bool noisyKsi1)
{
	clog<<"in calE2PCFsTrap in COSEBIs"<<endl;
	if(isFlat)
		omega_lambda=1.-omega_matter;
	if(nBinsTrapBegin<2)
		nBinsTrapBegin=2;

	NoiseKsi=0.;
	noisyKsi=noisyKsi1;
	if(noisyKsi)
	{
		clog<<"Adding noise to Ksi"<<endl;
		random.setParameters(0.,1.);
	}
	int nSteps= (nBinsTrapEnd-nBinsTrapBegin)/stepSize+2;
	//clog<<"#nSteps="<<nSteps<<endl;
	number deltaTheta=thetamax-thetamin;
	//cout<<"#deltaTheta="<<deltaTheta<<endl;
	matrix En(nMaximum*nPairs+1,nSteps);
	matrix EnP(nMaximum*nPairs+1,nSteps);
	matrix EnM(nMaximum*nPairs+1,nSteps);
	setTs(nMaximum);
	setKsi(sigma8,omega_matter,omega_lambda,w0,nindex,h100,omega_baryon);

	for(int i=0; i<nSteps; i++)
	{
		int nBinsTrap=nBinsTrapBegin+stepSize*i;
		clog<<"#nBinsTrap="<<nBinsTrap<<endl;
		number h=deltaTheta/(nBinsTrap-1.);
		//clog<<"#h="<<h<<endl;
		if(noisyKsi)
		{
			number dTheta=h;
			NoiseKsi=noise_vec[0]/sqrt(A*pi*dTheta/2);
			random_matP.resize(nBinsTrap);
			random_matM.resize(nBinsTrap);
			for(int nT=0;nT<nBinsTrap;nT++)
			{
				random_matP.load(nT,random.get());
				random_matM.load(nT,random.get());
			}
			//clog<<"NoiseKsi="<<NoiseKsi<<endl;
		}
		En.load(0,i,nBinsTrap);
		EnP.load(0,i,nBinsTrap);
		EnM.load(0,i,nBinsTrap);
		
		
		for(int bin1=0;bin1<nBins;bin1++)
		{
			for(int bin2=bin1;bin2<nBins;bin2++)
			{
				int m=0;
				int p1=calP(nBins,bin1,bin2);
				redshiftPair=p1;
				for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
				{
					//clog<<"loading E n="<<n1+1<<endl;
					matrix E3=valueEn2PCFsTrap(m,nBinsTrap,h);
					En.load(n1+1,i,E3.get(0));
					EnP.load(n1+1,i,E3.get(1));
					EnM.load(n1+1,i,E3.get(2));
				}
			}
		}
	}
	vector<matrix> En_vec;
	En_vec.push_back(En);
	En_vec.push_back(EnP);
	En_vec.push_back(EnM);
	return En_vec;
}*/

///the noise part only works for no tomography for now and is very simplified
/*matrix COSEBIs::calEn2PCFsTrap(int nBinsTrap,bool noisyKsi1)
{
	clog<<"in calE2PCFsTrap with a given nTrap in COSEBIs"<<endl;

	NoiseKsi=0.;
	noisyKsi=noisyKsi1;
	if(noisyKsi)
	{
		clog<<"Adding noise to Ksi"<<endl;
		random.setParameters(0.,1.);
	}
	number deltaTheta=thetamax-thetamin;
	matrix En(nMaximum*nPairs);
	setTs(nMaximum);
	setKsi(sigma8,omega_matter,omega_lambda,w0,nindex,h100,omega_baryon);
	
		

	clog<<"#nBinsTrap="<<nBinsTrap<<endl;
	number h=deltaTheta/(nBinsTrap-1.);
	if(noisyKsi)
	{
		number dTheta=h;
		NoiseKsi=noise_vec[0]/sqrt(A*pi*dTheta/2);
		clog<<"NoiseKsi="<<NoiseKsi<<endl;
		random_matP.resize(nBinsTrap);
		random_matM.resize(nBinsTrap);
		for(int nT=0;nT<nBinsTrap;nT++)
		{
			random_matP.load(nT,random.get());
			random_matM.load(nT,random.get());
		}
	}
	
	for(int bin1=0;bin1<nBins;bin1++)
	{
		for(int bin2=bin1;bin2<nBins;bin2++)
		{
			int m=0;
			int p1=calP(nBins,bin1,bin2);
			redshiftPair=p1;
			for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
			{
				matrix E3=valueEn2PCFsTrap(m,nBinsTrap,h);
				En.load(n1,E3.get(0));
			}
		}
	}
	return En;
}*/

// matrix COSEBIs::valueEn2PCFsKsiInput(int n,number h)
// {
// 	Cov_on=false;
// 	Real=true;
// 	number IntPlus= 0.;
// 	number IntMinus= 0.;
// 	nT=n-1;
// 	//cout<<"#n="<<n<<endl;
// 	realPlus=true;
// 	iRand=0;
// 	number integThetaminPlus=trapIntegrant(MinTheta);
// 	iRand=nBinsTrap-1;
// 	number integThetamaxPlus=trapIntegrant(MaxTheta);
// 	//cout<<"#integThetaminPlus="<<integThetaminPlus<<endl;
// 	//cout<<"#integThetamaxPlus="<<integThetamaxPlus<<endl;
// 	realPlus=false;
// 	iRand=0;
// 	number integThetaminMinus=trapIntegrant(thetamin);
// 	iRand=nBinsTrap-1;
// 	number integThetamaxMinus=trapIntegrant(thetamax);
// 	//cout<<"#integThetaminMinus="<<integThetaminMinus<<endl;
// 	//cout<<"#integThetamaxMinus="<<integThetamaxMinus<<endl;
// 	for(int i=1;i<(nBinsTrap-1);i++)
// 	{
// 		//clog<<"i="<<i<<endl;
// 		iRand=i;
// 		realPlus=true;
// 		number t=h*number(i)+thetamin;
// 		number IntP=integrant(t);
// 		realPlus=false;
// 		number IntM=integrant(t);
// 		IntMinus+=IntM;
// 		IntPlus+=IntP;
// 		//cout<<i<<"\t"<<t<<"\t"<<IntPlus<<"\t"<<IntMinus<<endl;
// 	}
// 
// 	IntPlus+=integThetaminPlus/2.;
// 	IntPlus+=integThetamaxPlus/2.;
// 	IntMinus+=integThetaminMinus/2.;
// 	IntMinus+=integThetamaxMinus/2.;
// 	IntPlus=IntPlus*h;
// 	IntMinus=IntMinus*h;
// 	matrix E3(3);
// 	//E3.load(0,n);
// 	E3.load(0,(IntPlus+IntMinus)/2.);
// 	E3.load(1,IntPlus);
// 	E3.load(2,IntMinus);
// 	return E3;
// }


// ///finds the value of xi_E/B+ for a given theta and d theta=h
// matrix COSEBIs::valuexiEB2PCFsKsiInput(matrix& Ksi_mat,number theta,number h)
// {
// 	number IntPrime1= 0.;// int_theta^inf d theta/theta  xi_-
// 	number IntPrime2= 0.;// int_theta^inf d theta/theta^3  xi_-
// 	int nTheta=Ksi_mat.rows;
// 	//ksi_mat[0]: theta, ksi_mat[1]: xi_+, ksi_mat[2]: xi_-
// 	//matrix intp_mat(2,nTheta);
// 	//matrix intm_mat(2,nTheta);
// 	number xi1,xi2;
// 	// 	Ksi_mat.printOut((string("Ksi_mat_2.ascii")).c_str(),10);
// 	// 	clog<<"in valueEn2PCFsKsiInput, nTheta="<<nTheta<<endl;
// 	intp=(Ksi_mat.get(0,0)*Ksi_mat.get(1,0)*Tpm_mat_vec[m].get(0,0))/2.;
// 	IntPlus=intp;
// 	//intp_mat.load(0,0,Ksi_mat.get(0,0));
// 	//intp_mat.load(1,0,intp);
// 	// 	clog<<"IntPlus="<<IntPlus<<", theta="<<Ksi_mat.get(0,0)<<", Ksip="<<Ksi_mat.get(1,0)
// 	// 		<<", Tp="<<Tpm_mat_vec[m].get(0,0)<<endl;
// 	intp=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(1,nTheta-1)*Tpm_mat_vec[m].get(0,nTheta-1))/2.;
// 	IntPlus+=intp;
// 	//intp_mat.load(0,nTheta-1,Ksi_mat.get(0,nTheta-1));
// 	//intp_mat.load(1,nTheta-1,intp);
// 	//clog<<"IntPlus="<<IntPlus<<", theta="<<Ksi_mat.get(0,nTheta-1)<<", Ksip="<<Ksi_mat.get(1,nTheta-1)
// 	//	<<", Tp="<<Tpm_mat_vec[m].get(0,nTheta-1)<<endl;
// 	intm=(Ksi_mat.get(0,0)*Ksi_mat.get(2,0)*Tpm_mat_vec[m].get(1,0))/2.;
// 	IntMinus=intm;
// 	//intm_mat.load(0,nTheta-1,Ksi_mat.get(0,nTheta-1));
// 	//intm_mat.load(1,nTheta-1,intm);
// 	IntMinus+=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(2,nTheta-1)*Tpm_mat_vec[m].get(1,nTheta-1))/2.;
// 	// 	clog<<IntPlus<<"\t"<<IntMinus<<endl;
// 	for(int itheta=1; itheta<(nTheta-1); itheta++)
// 	{
// 		xiprime1=Ksi_mat.get(0,itheta)
// 		intp=Ksi_mat.get(0,itheta)*Ksi_mat.get(1,itheta)*Tpm_mat_vec[m].get(0,itheta);
// 		IntPlus +=intp;
// 		intm=Ksi_mat.get(0,itheta)*Ksi_mat.get(2,itheta)*Tpm_mat_vec[m].get(1,itheta);
// 		IntMinus+=intm;
// 	}
// 
// 	IntPlus*=h;
// 	IntMinus*=h;
// 	number IntTotal=(IntPlus+IntMinus)/2.;
// 	
// 	matrix E3(3);
// 	E3.load(0,IntTotal);
// 	E3.load(1,IntPlus);
// 	E3.load(2,IntMinus);
// 	return E3;
// }

// matrix COSEBIs::calxiEB2PCFsFromInputKsi(vector<string> FileName,vector<string> corrFile,int nColumns)
// {
// 	clog<<"in calxiEB2PCFsFromInputKsi in COSEBIs with corrections"<<endl;
// 	for(unsigned int i=0; i<FileName.size(); i++)
// 	{
// 		clog<<"FileName is:"<<FileName[i]<<endl;
// 		clog<<"Correction FileName is:"<<corrFile[i]<<endl;
// 	}
// 	
// 	//number deltaTheta=thetamax-thetamin;
// // 	matrix En(5,nMaximum*nPairs);
// // 	setTs(nMaximum);
// 	//number of redshift pairs
// 	nPairs=FileName.size();
// 	//clog<<"nPairs="<<nPairs<<endl;
// 	nBins=(sqrt(8*nPairs+1)-1)/2;
// 	for(int r=0;r<nPairs;r++)
// 	{
// 		clog<<"r="<<r<<endl;
// 		vector<vector<number> > Ksi_vec=readKsi(FileName[r],nColumns);
// 		vector<vector<number> > Corr_vec=readKsi(corrFile[r],nColumns);
// 		// 		clog<<Corr_vec[0][0]<<endl;
// 		
// 		vector<int> index=FindMinMaxTheta(Ksi_vec);
// 		//clog<<"index[0]="<<index[0]<<"  index[1]="<<index[1]<<endl;
// 		matrix theta_mat(index[1]-index[0]+1);
// 		matrix Ksi_mat(3,index[1]-index[0]+1);
// 		//matrix Ksi_matnoCorr(3,index[1]-index[0]+1);
// 		//matrix Corr_mat(2,index[1]-index[0]+1);
// 		for(int i=index[0];i<=index[1]; i++)
// 		{
// 			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
// 			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
// 			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[i][1]);
// 			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[i][1]);
// 		}
// 		// 		Ksi_mat.printOut((folderName+string("/Ksi_mat.ascii")).c_str(),10);
// 		// 		exit(1);
// 		//Ksi_matnoCorr.printOut((folderName+string("/Ksi_matnoCorr.ascii")).c_str(),10);
// 		//Corr_mat.printOut((folderName+string("/Corr_mat.ascii")).c_str(),10);
// 		
// 		//exit(1);		
// 		//clog<<"theta1="<<Ksi_mat.get(0,0)/arcmin<<"  theta2="<<Ksi_mat.get(0,1)/arcmin<<endl;
// 		//clog<<"Ksip(theta1)="<<Ksi_mat.get(1,0)<<"   Ksip(theta2)="<<Ksi_mat.get(1,1)<<endl;
// 		//clog<<"Ksim(theta1)="<<Ksi_mat.get(2,0)<<"   Ksim(theta2)="<<Ksi_mat.get(2,1)<<endl;
// 		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
// 		//clog<<"h="<<h<<endl;
// // 		FindTnKsiInput(theta_mat);
// 		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
// 		{
// 			matrix xi=valueEn2PCFsKsiInput(Ksi_mat,m,h);
// 			xiEB.load(0,n1,theta);
// 			xiEB.load(1,n1,xi.get(0));
// 			xiEB.load(2,n1,xi.get(1));
// 			xiEB.load(3,n1,xi.get(2));
// 		}
// 	}
// 	return xiEB;
// }






/*
matrix COSEBIs::Derivative(matrix& EnM,number h,int size)
{
	matrix  derivativeE(size);
	number der;
	for (int n=0;n<size;n++)
	{
		if(derivative==2)
			der=(-EnM.get(0,n)+EnM.get(1,n))/(2.*h);
		else if(derivative==4)
			der=(EnM.get(0,n)-8.*EnM.get(1,n)+8.*EnM.get(2,n)-EnM.get(3,n))/(12.*h);
		else if(derivative==6)
			der=(-EnM.get(0,n)+9.*EnM.get(1,n)-45.*EnM.get(2,n)
				+45.*EnM.get(3,n)-9.*EnM.get(4,n)+EnM.get(5,n))/(60.*h);
		derivativeE.load(n,der);
	}
	return derivativeE;
}

void COSEBIs::setParamValueDEn(int i,int par,number h)
{
  clog<<"COSEBIs:in setParamValueDEn: setting the parameter's value derivative is:"
		<<derivative<<endl;
	if(derivative==2)
		switch(i)
		{
			case 0:
				param_vec[par]-=h;
			break;
			case 1:
				param_vec[par]+=2*h;
			break;
		}
	else if(derivative==4)
		switch(i)
		{
			case 0:
				param_vec[par]-=2.*h;
			break;
			case 1:
				param_vec[par]+=h;
			break;
			case 2:
				param_vec[par]+=2.*h;
			break;
			case 3:
				param_vec[par]+=h;
			break;
		}
	else if (derivative==6)
		switch(i)
		{
			case 0:
				param_vec[par]-=3.*h;
			break;
			case 1:
				param_vec[par]+=h;
			break;
			case 2:
				param_vec[par]+=h;
			break;
			case 3:
				param_vec[par]+=2.*h;
			break;
			case 4:
				param_vec[par]+=h;
			break;
			case 5:
				param_vec[par]+=h;
			break;
		}
	clog.precision(6);
	clog<<"prameter "<<par<<" value is:"<<param_vec[par]<<endl;
}

matrix COSEBIs::calDEn()
{
	clog<<"in calDEn in COSEBIs"<<endl;
	setWns(nMaximum);
	determine_integration_limits_En();
	DEn_mat.resize(param, nPairs*nMaximum);
	string FlatString="";
	if(isFlat)
		FlatString="_Flat";
// 		omega_lambda=1.-omega_matter;
	string paramstring="_";
	for(int i=0; i<param; i++)
	{
		paramstring+=paramsOrder[i]+"_";
	}
	string DEnFileName=folderName+string("/DEn_max-")
				+powerswitchName+FlatString+string("-")
				+toString(derivative)+string("_nPar_")+toString(param)+paramstring
				+string("nMax_")+toString(nMaximum)+string("_Z_")
				+toString(begin,2)+string("-")+toString(end,2)+string("_nBins_")
				+toString(nBins)+string("-s8_")
				+toString(sigma8,4)+string("-Om_")
				+toString(omega_matter,4)+string("-OL_")
				+toString(omega_lambda,4)+string("-W0_")
				+toString(w0,4)+string("-nindex_")
				+toString(nindex,4)+string("-h100_")
				+toString(h100,4)+string("-Ob_")
				+toString(omega_baryon,4)+string("_th_")
				+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)+string(".bin");
	ifstream DEnfile(DEnFileName.c_str());
	if(DEnfile.fail())
	{
		clog<<"DEn Failed:"<<DEnFileName<<endl;
		initializeParam();
		matrix Eparam_mat(derivative,nMaximum*nPairs);
		vector<matrix> DEparam_vec;
		for(int i=0; i<param; i++)
			DEparam_vec.push_back(matrix(nMaximum*nPairs));
		vector<vector<number> >z_sample(2, vector<number> (2));
		//clog<<"DEn vec made"<<endl;
		//the number of the parameter considered
		for(int par=0;par<param;par++)
		{
				//clog<<"in par loop par="<<par<<endl;
			number h=param_vec[par]*H;
			clog<<"parameter("<<par<<")="<<param_vec[par]<<endl;
			//initializing the 4 values for each set for derivation
			for(int i=0; i<derivative;i++)
			{
				//sets the value of the param_vec for each i to 
				//initialize the En vector for four point derivative
				setParamValueDEn(i,par,h);
				calPower(param_vec, paramsOrder);
				for(int bin1=0;bin1<nBins;bin1++)
				{
					for(int bin2=bin1;bin2<nBins;bin2++)
					{
						int m=0;
						int p1=calP(nBins,bin1,bin2);
						redshiftPair=p1;
						for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
						{
							clog<<"in DEn loading E n="<<n1+1<<endl;
							Eparam_mat.load(i,n1,valueEn(m));
						}
					}
				}
				if(i==(derivative-1))
					initializeParam();
			}//end of for(i)
			DEparam_vec[par]=Derivative(Eparam_mat,h,nMaximum*nPairs);
		}
		for(int p=0;p<param;p++)
			for(int n=0; n<nPairs*nMaximum; n++)
				DEn_mat.load(p,n,DEparam_vec[p].get(n));
		DEn_mat.store(DEnFileName.c_str());
	}
	else
		DEn_mat.restore(DEnFileName.c_str());
	DEnfile.close();
	DEn_calculated=true;
	return DEn_mat;
}

matrix COSEBIs::calM(int i, int j, int nMax)
{
	if(!DEn_calculated)
		DEn_mat=calDEn();
	//clog<<"in COSEBIs in calM"<<endl;
	matrix M(nMax*nPairs,nMax*nPairs);
	matrix DEParami(nMax*nPairs);
	matrix DEParamj(nMax*nPairs);
	for(int n=0; n<nMax*nPairs; n++)
	{
		int m=ceil(n/nMax)*nMaximum+n%nMax;
		//clog<<"m="<<m<<endl;
		DEParami.load(n,DEn_mat.get(i,m));
		DEParamj.load(n,DEn_mat.get(j,m));
	}
	//clog<<"DEparami and DEparamj are done"<<endl;
	M=DEParami*DEParamj.t()+DEParamj*DEParami.t();
  	return M;
}
///only for flat cosmology
matrix COSEBIs::calMW0Last(int i, int j, int nMax)
{
	if(!DEn_calculated)
		DEn_mat=calDEn();
	//clog<<"in COSEBIs in calM"<<endl;
	matrix M(nMax*nPairs,nMax*nPairs);
	matrix DEParami(nMax*nPairs);
	matrix DEParamj(nMax*nPairs);
	for(int n=0; n<nMax*nPairs; n++)
	{
		int m=ceil(n/nMax)*nMaximum+n%nMax;
		//clog<<"m="<<m<<endl;
		if(i==0 || i==1)
			DEParami.load(n,DEn_mat.get(i,m));
		else if(i<5)
			DEParami.load(n,DEn_mat.get(i+1,m));
		else
			DEParami.load(n,DEn_mat.get(2,m));
		

		if(j==0 || j==1)
			DEParamj.load(n,DEn_mat.get(j,m));
		else if(j<5)
			DEParamj.load(n,DEn_mat.get(j+1,m));
		else
			DEParamj.load(n,DEn_mat.get(2,m));
		
	}
	//clog<<"DEparami and DEparamj are done"<<endl;
	M=DEParami*DEParamj.t()+DEParamj*DEParami.t();
  	return M;
}

void COSEBIs::initializeParam()
{
	clog<<"inside initializeParam in COSEBIs"<<endl;
	param_vec.clear();
	for (int i=0; i<param;i++)
	{
		param_vec.push_back(0.);
		if (paramsOrder[i]=="sigma8")
			param_vec[i]=sigma8;
		else if (paramsOrder[i]=="omega_matter")
			param_vec[i]=omega_matter;
		else if (paramsOrder[i]=="omega_lambda")
			param_vec[i]=omega_lambda;
		else if (paramsOrder[i]=="nindex")
			param_vec[i]=nindex;
		else if (paramsOrder[i]=="w0")
			param_vec[i]=w0;
		else if (paramsOrder[i]=="h100")
			param_vec[i]=h100;
		else if (paramsOrder[i]=="omega_baryon")
			param_vec[i]=omega_baryon;
		else
		{
			clog<<paramsOrder[i]<<" is not a known parameter"<<endl;
			exit(1);
		}
	}
}

void COSEBIs::setParamValueDDEn(int i, int par1,int par2,number h,number k)
{
	clog<<"in set param value for second derivative COSEBIs"<<endl;
	clog.precision(5);
	//clog<<"h="<<h<<endl;
	if(par1-par2)
	{
		switch(i)
		{
			case 0:
				param_vec[par1]-=2.*h;
				param_vec[par2]-=2.*k;
			break;
			case 1:
				param_vec[par2]+=k;
			break;
			case 2:
				param_vec[par2]+=2.*k;
			break;
			case 3:
				param_vec[par2]+=k;
			break;

			case 4:
				param_vec[par1]+=h;
				param_vec[par2]-=4.*k;
			break;
			case 5:
				param_vec[par2]+=k;
			break;
			case 6:
				param_vec[par2]+=2.*k;
			break;
			case 7:
				param_vec[par2]+=k;
			break;

			case 8:
				param_vec[par1]+=2.*h;
				param_vec[par2]-=4.*k;
			break;
			case 9:
				param_vec[par2]+=k;
			break;
			case 10:
				param_vec[par2]+=2.*k;
			break;
			case 11:
				param_vec[par2]+=k;
			break;

			case 12:
				param_vec[par1]+=h;
				param_vec[par2]-=4.*k;
			break;
			case 13:
				param_vec[par2]+=k;
			break;
			case 14:
				param_vec[par2]+=2.*k;
			break;
			case 15:
				param_vec[par2]+=k;
			break;
		}
		clog<<"i="<<i<<endl;
		clog<<"h="<<h<<"   k="<<k<<endl;
		clog<<"param_vec["<<par1<<"]="<<param_vec[par1]<<"  param_vec["<<par2<<"]="<<param_vec[par2]<<endl;
	}
	else
	{
		switch(i)
		{
			case 0:
				param_vec[par1]-=2.*h;
			break;
			case 1:
				param_vec[par1]+=h;
			break;
			case 2:
				param_vec[par1]+=h;
			break;
			case 3:
				param_vec[par1]+=h;
			break;
			case 4:
				param_vec[par1]+=h;
		}
	clog<<"i="<<i<<endl;
	clog<<"h="<<h<<endl;
	clog<<"param_vec["<<par1<<"]="<<param_vec[par1]<<endl;
	}
}


void COSEBIs::setParamValueDDEnSimple(int i, int par1,int par2,number h,number k)
{
	clog<<"in set param value for simple second derivative COSEBIs"<<endl;
	clog.precision(5);
	//clog<<"h="<<h<<endl;
	if(par1-par2)
	{
		switch(i)
		{
			case 0:
				param_vec[par1]-=h;
				param_vec[par2]-=k;
			break;
			case 1:
				param_vec[par2]+=k;
			break;
			case 2:
				param_vec[par1]+=h;
				param_vec[par2]-=k;
			break;
			case 3:
				param_vec[par2]+=k;
			break;
			case 4:
				param_vec[par2]+=k;
			break;
			case 5:
				param_vec[par1]+=h;
				param_vec[par2]-=k;
			break;
			case 6:
				param_vec[par2]+=k;
			break;
		}
		clog<<"i="<<i<<endl;
		clog<<"h="<<h<<"   k="<<k<<endl;
		clog<<"param_vec["<<par1<<"]="<<param_vec[par1]<<"  param_vec["<<par2<<"]="<<param_vec[par2]<<endl;
	}
	else
	{
		switch(i)
		{
			case 0:
				param_vec[par1]+=h;
			break;
			case 1:
				param_vec[par1]-=h;
			break;
			case 2:
				param_vec[par1]-=h;
			break;
		}
	clog<<"i="<<i<<endl;
	clog<<"h="<<h<<endl;
	clog<<"param_vec["<<par1<<"]="<<param_vec[par1]<<endl;
	}
}

matrix COSEBIs::SecondDerivative(matrix & En_mat,number h, number k)
{
	matrix DD_mat(nMaximum*nPairs);
	if(OneParam)
	{
		for(int n=0; n<nMaximum*nPairs; n++)
		{
		number result=(En_mat.get(0,n)-2*En_mat.get(1,n)+En_mat.get(2,n))/(h*h);
		DD_mat.load(n,result);
		}	
	}
	else
		for(int n=0; n<nMaximum*nPairs; n++)
		{
			number result=(En_mat.get(0,n)-En_mat.get(1,n)-En_mat.get(2,n)+2*En_mat.get(3,n)
								-En_mat.get(4,n)-En_mat.get(5,n)+En_mat.get(6,n))/(2.*k*h);
			DD_mat.load(n,result);
		}
	return DD_mat;
}


matrix COSEBIs::fourPointSecondDerivative(matrix& En_mat,number h,number k)
{
	matrix DD_mat(nMaximum*nPairs);
	if(OneParam)
	{
		for(int n=0; n<nMaximum*nPairs; n++)
		{
			number result=(-En_mat.get(0,n)
								+16.*En_mat.get(1,n)
								-30.*En_mat.get(2,n)
								+16.*En_mat.get(3,n)
								-En_mat.get(4,n))/(12.*h*h);
			DD_mat.load(n,result);
		}
	}
	else
	{
		for(int n=0; n<nMaximum*nPairs; n++)
		{
			number	 result=1./(600.*h*k)*
					(-63.*(En_mat.get(8,n)+En_mat.get(13,n)+En_mat.get(2,n)+En_mat.get(7,n))
					+63.*(En_mat.get(4,n)+En_mat.get(1,n)+En_mat.get(11,n)+En_mat.get(14,n))
					+44.*(En_mat.get(12,n)+En_mat.get(3,n)-En_mat.get(0,n)-En_mat.get(15,n))
					+74.*(En_mat.get(5,n)+En_mat.get(10,n)-En_mat.get(9,n)-En_mat.get(6,n)));
			DD_mat.load(n,result);
		}
	}
	return DD_mat;
}

vector<vector<matrix> > COSEBIs::calDDE()
{
	clog<<"calculating DDEn in COSEBIs"<<endl;
	Cov_on=false;
	setWns(nMaximum);
	determine_integration_limits_En();
	vector< vector< number> > z_sample(2, vector<number> (2));
	matrix En_mat1(5,nMaximum*nPairs);
	matrix En_mat2(16,nMaximum*nPairs);
	vector<number> h_vec(param);
	initializeParam();
	for(int i=0; i<param; i++)
	{
		vector<matrix> DDEn_vec;
		for(int j=0; j<param; j++)
			DDEn_vec.push_back(matrix(nMaximum*nPairs));
		DDEn_vecvec.push_back(DDEn_vec);
		h_vec[i]=H*param_vec[i];
	}
	OneParam=true;
	string FlatString="";
	if(isFlat)
		FlatString="_Flat";
	for(int par1=0; par1<param; par1++)
	{
	  	string DDEnFileName=folderName+string("/DDEn-")
				+powerswitchName+FlatString+string("-p")
				+toString(par1)+string("-p")+toString(par1)
				+string("_")+toString(nMaximum)+string("_")
				+toString(begin,2)+string("-")+toString(end,2)+string("_")
				+toString(nBins)+string("-")
				+toString(sigma8,4)+string("-")
				+toString(omega_matter,4)+string("-")
				+toString(omega_lambda,4)+string("-")
				+toString(w0,4)+string("-")
				+toString(nindex,4)+string("-")
				+toString(h100,4)+string("-")
				+toString(omega_baryon,4)+string("_")
				+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)+string(".bin");
		ifstream DDEnfile((DDEnFileName).c_str());
		if(DDEnfile.fail())
		{
		  	clog<<"DDEnfile Failed:"<<DDEnFileName<<endl;
			initializeParam();
			if(isFlat)
				for(int d=0; d<5; d++)
				{
					setParamValueDDEn(d,par1,par1,h_vec[par1],h_vec[par1]);
// 					AssessParamOrderAndChange(param_vec, paramsOrder);
					calPower(param_vec, paramsOrder);
					for(int bin1=0;bin1<nBins;bin1++)
					{
						for(int bin2=bin1;bin2<nBins;bin2++)
						{
							int m=0;
							int p1=calP(nBins,bin1,bin2);
							redshiftPair=p1;
							for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
							{
								clog<<"in DDEn loading E n="<<n1+1<<endl;
								En_mat1.load(d,n1,valueEn(m));
							}
						}
					}
				}
			DDEn_vecvec[par1][par1]=fourPointSecondDerivative(En_mat1,h_vec[par1],h_vec[par1]);
			DDEn_vecvec[par1][par1].store((DDEnFileName).c_str());
		}
		else
			DDEn_vecvec[par1][par1].restore((DDEnFileName).c_str());
		DDEnfile.close();
// 		DDEn_vecvec[par1][par1].printOut((folderName+string("/DDEn-")
// 				+powerswitchName+string("-p")
// 				+toString(par1)+string("-p")+toString(par1)
// 				+string("_")+toString(nMaximum)+string(".ascii")).c_str(),8);
	}
	if(param>1)
	{
		clog<<"in param>1 DDEn"<<endl;
		OneParam=false;
		for(int par1=0; par1<param; par1++)
		{
			for(int par2=par1+1; par2<param; par2++)
			{
			  	string DDEnFileName=folderName+string("/DDEn-")
							+powerswitchName+FlatString+string("-p")
							+toString(par1)
							+string("-p")+toString(par2)+string("_")
							+toString(nMaximum)+string("_")
							+toString(begin,2)+string("-")+toString(end,2)+string("_")
							+toString(nBins)+string("-")
							+toString(sigma8,4)+string("-")
							+toString(omega_matter,4)+string("-")
							+toString(omega_lambda,4)+string("-")
							+toString(w0,4)+string("-")
							+toString(nindex,4)+string("-")
							+toString(h100,4)+string("-")
							+toString(omega_baryon,4)+string("_")
							+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)
							+string(".bin");
				ifstream DDEnfile((DDEnFileName).c_str());
				if(DDEnfile.fail())
				{
					clog<<"DDEnfile Failed:"<<DDEnFileName<<endl;
					initializeParam();
					for(int d=0; d<16; d++)
					{
							setParamValueDDEn(d,par1,par2,h_vec[par1],h_vec[par2]);
// 							AssessParamOrderAndChange(param_vec, paramsOrder);
							calPower(param_vec, paramsOrder);
							for(int bin1=0;bin1<nBins;bin1++)
							{
								for(int bin2=bin1;bin2<nBins;bin2++)
								{
									int m=0;
									int p1=calP(nBins,bin1,bin2);
									redshiftPair=p1;
									for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
									{
										clog<<"in DDEn param>1 loading E n="<<n1+1<<endl;
										En_mat2.load(d,n1,valueEn(m));
									}
								}
							}
					}
					DDEn_vecvec[par1][par2]=fourPointSecondDerivative(En_mat2,h_vec[par1],h_vec[par2]);
					DDEn_vecvec[par1][par2].store((DDEnFileName).c_str());
				}
				else
					DDEn_vecvec[par1][par2].restore((DDEnFileName).c_str());

				DDEn_vecvec[par2][par1]=DDEn_vecvec[par1][par2];
				DDEnfile.close();
// 				DDEn_vecvec[par1][par2].printOut((folderName+string("/DDEn-")
// 						+powerswitchName+string("-p")
// 						+toString(par1)+string("-p")+toString(par2)
// 						+string("_")+toString(nMaximum)+string(".ascii")).c_str(),8);
			}
		}
	}
	return DDEn_vecvec;
}

vector<vector<matrix> > COSEBIs::calDDESimple()
{
	clog<<"calculating DDEn Simple in COSEBIs"<<endl;
	Cov_on=false;
	setWns(nMaximum);
	determine_integration_limits_En();
	vector< vector< number> > z_sample(2, vector<number> (2));
	matrix En_mat1(3,nMaximum*nPairs);
	matrix En_mat2(7,nMaximum*nPairs);
	vector<number> h_vec(param);
	initializeParam();
	for(int i=0; i<param; i++)
	{
		vector<matrix> DDEn_vec;
		for(int j=0; j<param; j++)
			DDEn_vec.push_back(matrix(nMaximum*nPairs));

		DDEn_vecvec.push_back(DDEn_vec);
		h_vec[i]=H*param_vec[i];
	}
	OneParam=true;
	string FlatString="";
	if(isFlat)
		FlatString="_Flat";
	for(int par1=0; par1<param; par1++)
	{
	  	string DDEnFileName=folderName+string("/DDEnSimple-")
				+powerswitchName+FlatString+string("_")+paramsOrder[par1]+string("-")+paramsOrder[par1]
				+string("_")+toString(nMaximum)+string("_")
				+toString(begin,2)+string("-")+toString(end,2)+string("_")
				+toString(nBins)+string("-")
				+toString(sigma8,4)+string("-")
				+toString(omega_matter,4)+string("-")
				+toString(omega_lambda,4)+string("-")
				+toString(w0,4)+string("-")
				+toString(nindex,4)+string("-")
				+toString(h100,4)+string("-")
				+toString(omega_baryon,4)+string("_")
				+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)+string(".bin");
		ifstream DDEnfile((DDEnFileName).c_str());
		if(DDEnfile.fail())
		{
		  	clog<<"DDEnfile Failed:"<<DDEnFileName<<endl;
			initializeParam();
			for(int d=0; d<3; d++)
			{
				setParamValueDDEnSimple(d,par1,par1,h_vec[par1],h_vec[par1]);
// 				AssessParamOrderAndChange(param_vec, paramsOrder);
				calPower(param_vec, paramsOrder);
				for(int bin1=0;bin1<nBins;bin1++)
				{
					for(int bin2=bin1;bin2<nBins;bin2++)
					{
						int m=0;
						int p1=calP(nBins,bin1,bin2);
						redshiftPair=p1;
						for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
						{
							clog<<"in DDEn loading E n="<<n1+1<<endl;
							En_mat1.load(d,n1,valueEn(m));
						}
					}
				}
			}
// 			En_mat1.store((folderName+string("/En_mat1-")+toString(par1)
// 				+string("-")+toString(nMaximum)+string("-")
// 				+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)
// 				+string(".ascii")).c_str(),10);
			DDEn_vecvec[par1][par1]=SecondDerivative(En_mat1,h_vec[par1],h_vec[par1]);
			DDEn_vecvec[par1][par1].store((DDEnFileName).c_str());
		}
		else
			DDEn_vecvec[par1][par1].restore((DDEnFileName).c_str());
		DDEnfile.close();
// 		DDEn_vecvec[par1][par1].printOut((folderName+string("/DDEnSimple-")
// 				+powerswitchName+string("-p")
// 				+toString(par1)+string("-p")+toString(par1)
// 				+string("_")+toString(nMaximum)+string(".ascii")).c_str(),8);
	}

	if(param>1)
	{
		clog<<"in param>1 DDEn"<<endl;
		OneParam=false;
		for(int par1=0; par1<param; par1++)
		{
			for(int par2=par1+1; par2<param; par2++)
			{
			  	string DDEnFileName=folderName+string("/DDEnSimple-")
							+powerswitchName+FlatString+string("_")+paramsOrder[par1]
							+string("-")+paramsOrder[par2]+string("_")
							+toString(nMaximum)+string("_")
							+toString(begin,2)+string("-")+toString(end,2)+string("_")
							+toString(nBins)+string("-")
							+toString(sigma8,4)+string("-")
							+toString(omega_matter,4)+string("-")
							+toString(omega_lambda,4)+string("-")
							+toString(w0,4)+string("-")
							+toString(nindex,4)+string("-")
							+toString(h100,4)+string("-")
							+toString(omega_baryon,4)+string("_")
							+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)+string(".bin");
				ifstream DDEnfile((DDEnFileName).c_str());
				if(DDEnfile.fail())
				{
					clog<<"DDEnfile Failed:"<<DDEnFileName<<endl;
					initializeParam();
					for(int d=0; d<7; d++)
						{
							setParamValueDDEnSimple(d,par1,par2,h_vec[par1],h_vec[par2]);
// 							AssessParamOrderAndChange(param_vec, paramsOrder);
							calPower(param_vec, paramsOrder);
							for(int bin1=0;bin1<nBins;bin1++)
							{
								for(int bin2=bin1;bin2<nBins;bin2++)
								{
									int m=0;
									int p1=calP(nBins,bin1,bin2);
									redshiftPair=p1;
									for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
									{
										clog<<"in DDEn param>1 loading E n="<<n1+1<<endl;
										En_mat2.load(d,n1,valueEn(m));
									}
								}
							}
						}
					DDEn_vecvec[par1][par2]=SecondDerivative(En_mat2,h_vec[par1],h_vec[par2]);
					DDEn_vecvec[par1][par2].store((DDEnFileName).c_str());
				}
				else
				{
					DDEn_vecvec[par1][par2].restore((DDEnFileName).c_str());
				}
				DDEn_vecvec[par2][par1]=DDEn_vecvec[par1][par2];
				DDEnfile.close();
// 				DDEn_vecvec[par1][par2].printOut((folderName+string("/DDEnSimple-")
// 						+powerswitchName+string("-p")
// 						+toString(par1)+string("-p")+toString(par2)
// 						+string("_")+toString(nMaximum)+string(".ascii")).c_str(),8);
			}
		}
	}
	return DDEn_vecvec;
}



matrix COSEBIs::calQ(matrix deltaPhi_est,int nMax)
{
	clog<<"calculating Q in COSEBIs"<<endl;
	matrix Q(param,nMax*nPairs);
	Q.zero();
	//cout<<"in calQ in COSEBIs"<<endl;
	for(int n=0; n<nMax*nPairs; n++)
		for(int p1=0; p1<param; p1++)
			for(int p2=0; p2<param; p2++)
			{
 				number result= DDEn_vecvec[p1][p2].get(ceil(n/nMax)*nMaximum+n%nMax)
									*deltaPhi_est.get(p2);
				Q.load(p1,n,(Q.get(p1,n)+result));
			}
	return Q;
}*/
  /*
vector<matrix> COSEBIs::returnMapStatisticsIntegrant(vector <number> ell_vec, number theta, int redshift)
{
	vector<matrix> IntegAndU_mat;
	matrix integ_mat(2, ell_vec.size());
	matrix U_mat(2, ell_vec.size());
	for(int i=0; i<ell_vec.size();i++)
	{
		number thetaell=theta/2.*ell_vec[i];
		number U_t=24.*gsl_sf_bessel_Jnu(4,thetaell)/(thetaell*thetaell);
		U_mat.load(0,i,ell_vec[i]);
		U_mat.load(1,i,U_t);
	}
	//calPower();
// 	nW=n-1;
	redshiftPair=redshift;
	for(int i=0; i<ell_vec.size();i++)
	{
		number integ=U_mat.get(1,i)*U_mat.get(1,i)*ell_vec[i]*powerspectrum_vec[redshiftPair].value(ell_vec[i]);
		integ_mat.load(0,i,ell_vec[i]);
		integ_mat.load(1,i,integ);
	}
	IntegAndU_mat.push_back(integ_mat);
	IntegAndU_mat.push_back(U_mat);
	return IntegAndU_mat;
}*/
