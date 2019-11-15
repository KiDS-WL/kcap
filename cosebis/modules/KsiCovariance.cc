#include "KsiCovariance.h"

KsiCovariance::KsiCovariance(){}
KsiCovariance::~KsiCovariance(){}
KsiCovariance::KsiCovariance(number thetamin,number thetamax,int nKsiBins)
{
	setParameters(thetamin,thetamax,nKsiBins);
}

// if(logbinning)
// {
// 	number theta_binned=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+0.5));
// 	// 							clog<<"logbining, theta_binned="<<theta_binned/arcmin<<endl;
// 	theta_max=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+1.));
// 	theta_min=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta));

// }
// else
// {
// 	number theta_binned=(tmax-tmin)/(nTheta_binned)*(itheta+0.5)+tmin;
// 	// 							clog<<"linbinning, theta_binned="<<theta_binned/arcmin<<endl;
// 	theta_max=theta_binned+(tmax-tmin)/(nTheta_binned)*0.5;
// 	theta_min=theta_binned-(tmax-tmin)/(nTheta_binned)*0.5;
// }


void KsiCovariance::setParameters(number thetamin1,number thetamax1,int nKsiBins1)
{
	MixedTerm=false;
	nKsiBins=nKsiBins1;
	thetamin=thetamin1;
	thetamax=thetamax1;
	LogCoef=log(thetamax/thetamin)/(nKsiBins);
}

void KsiCovariance::setParameters(matrix theta_plus,matrix theta_minus)
{
	MixedTerm=false;
	//nKsiBins=nKsiBins1;
	//thetamin=thetamin1;
	//thetamax=thetamax1;
	nKsiBins=theta_plus.rows;
	LogCoef=log(theta_plus.get(theta_plus.size()-1)/theta_plus.get(0))/(nKsiBins);
}


void KsiCovariance::setZbins(int nPairs1)
{
	nPairs=nPairs1;
	nBins=int((-1+sqrt(1+8*nPairs))/2);
	//clog<<"nPairs="<<nPairs<<" nBins="<<nBins<<endl;
}



void KsiCovariance::setNoise(number A1,vector<number> sigma_e_vec1,vector<number> nBar_vec)
{
	clog<<"setting noise in KsiCovariance"<<endl;
	A=A1;
	noise_vec.clear();
	sigma_e_vec=sigma_e_vec1;
	for(int bin=0; bin<sigma_e_vec.size(); bin++)
	{
		noise_vec.push_back(sigma_e_vec[bin]*sigma_e_vec[bin]/(2.*nBar_vec[bin]));
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
	}
}

void KsiCovariance::SetNoNoise()
{
	NoNoise=true;
}

void KsiCovariance::setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM)
{
	clog<<"setting Ksi in KsiCovarinace"<<endl;
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


void KsiCovariance::setPower(vector<number> log_ell,vector<vector<number> > InputPower)
{
	//clog<<"in setPower in COSEBIs.cc"<<endl;
	//check if InputPower.size()== nPairs
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
// 		string powerFileName=string("Power_")+toString(r);
//   		powerspectrum_vec[r].setName(powerFileName.c_str(),function_cosebis::NONAMECOUNTER);
//  		powerspectrum_vec[r].makeTable(0.15,100000.,200,true);
//  		powerspectrum_vec[r].saveTable();
	}
	MinPowerCov=exp(log_ell[0]);
	MaxPowerCov=exp(log_ell[log_ell.size()-1]);
}


void KsiCovariance::set_neff(vector<number> neff_vec1)
{
	clog<<"setting neff"<<endl;
	neff_vec=neff_vec1;
}


number KsiCovariance::ReturnPower(number ell,int rPair)
{
	return powerspectrum_vec[rPair].value(ell);
}


// matrix KsiCovariance::returnIntegrand(vector <number> ell_vec, number theta, int redshift)
// {
// 	matrix integ_mat(2, ell_vec.size());
	
// 	setWns(nMaximum);
// 	nW=n-1;
// 	redshiftPair=redshift;
// 	for(int i=0; i<ell_vec.size();i++)
// 	{
// 		integ_mat.load(0,i,ell_vec[i]);
// 		integ_mat.load(1,i,integrant(ell_vec[i]));
// 	}
// 	return integ_mat;
// }


number KsiCovariance::delta(int i1, int j1)
{
	return i1==j1? 1.: 0.;
}

number KsiCovariance::integrant(number l)
{
	number integ=0.;
	number mixedPower=0.;

	if(MixedTerm)
		mixedPower=powerspectrum_vec[rp1].value(l)*delta2noise
				+powerspectrum_vec[rp2].value(l)*delta1noise
				+powerspectrum_vec[rp3].value(l)*delta4noise
				+powerspectrum_vec[rp4].value(l)*delta3noise;
	else
		mixedPower=powerspectrum_vec[rp1].value(l)*powerspectrum_vec[rp2].value(l)
							+powerspectrum_vec[rp1].value(l)*delta2noise
							+powerspectrum_vec[rp2].value(l)*delta1noise
							+powerspectrum_vec[rp3].value(l)*powerspectrum_vec[rp4].value(l)
							+powerspectrum_vec[rp3].value(l)*delta4noise
							+powerspectrum_vec[rp4].value(l)*delta3noise;

	//integrands of covariance for bessel functions of order 0 or 4
// 	if((theta1==1.*arcmin)&&(theta2==1.*arcmin)&& (orderB1==0)&&(orderB2==0))
// 			cout<<l<<"\t"<<mixedPower<<endl;
	if((orderB1==0)&&(orderB2==0))
	{
		integ=l*gsl_sf_bessel_J0(l*theta1)*gsl_sf_bessel_J0(l*theta2)*mixedPower;
	}

	else if((orderB1==0)&&(orderB2==4))
	{
		integ= l*gsl_sf_bessel_J0(l*theta1)*gsl_sf_bessel_Jn(4,l*theta2)*mixedPower;
	}

	else if((orderB1==4)&&(orderB2==0))
	{
		integ= l*gsl_sf_bessel_Jn(4,l*theta1)*gsl_sf_bessel_J0(l*theta2)*mixedPower;
	}

	else if((orderB1==4)&&(orderB2==4))
	{
		integ= l*gsl_sf_bessel_Jn(4,l*theta1)*gsl_sf_bessel_Jn(4,l*theta2)*mixedPower;
	}
	else
	{
		clog<<"orderB1 or orderB2 not valid"<<endl;
		exit(1);
	}
	return integ;
}

///value of covariance
number KsiCovariance::value(number theta11,number theta21,int orderB11,int orderB21)
{
	int accuracyG=60;
	number lLim=MinPowerCov;
	number uLim=MaxPowerCov;
	//covariance=true;
	number resultG = 0.;
	number resultR=0.;
	//difference=0.;
	//bessel function order
	orderB1=orderB11;
	orderB2=orderB21;
	theta1=theta11;
	theta2=theta21;
	number arg=0.;
	number deltaTheta=0.;
	const number epsPP=1e-8;

	bool terminate=false;
	//number noise=0.;
// 	if(theta1==theta2)
// 		deltaTheta=theta1*(exp(log(thetamax/thetamin)/(nKsiBins-1))-1);

	
	//cout<<"#orderB1="<<orderB1<<"   orderB2="<<orderB2<<endl;
// 	if(theta1==theta2)
// 	{
// 		//?
// 		if(Logscale)
// 			deltaTheta=theta1*(exp(0.5*LogCoef)-exp(-0.5*LogCoef));
// 		else
// 			deltaTheta=(thetamax-thetamin)/(nKsiBins-1.0);
// 		//cout<<"# deltaTheta="<<deltaTheta<<endl;
// 	}

// 	number noiseTheta;
// 	if(NoNoise)
// 		noiseTheta=0.;
// 	else
// 		noiseTheta=2.*deltaNoise*noise*noise/(theta1*deltaTheta);
	//cout<<"# deltaTheta="<<deltaTheta<<endl;

	if((orderB1==0)&&(orderB2==0))
	{
		int iL=0;
		int iS=0;
		number tLarge=0.;
		number tSmall=0.;
		if(theta1>theta2)
		{
			tLarge=theta1;
			tSmall=theta2;
		}
		else
		{
			tLarge=theta2;
			tSmall=theta1;
		}
		//cout<<"#tLarge="<<tLarge<<"  tSmall="<<tSmall<<endl;
		//cout<<"# case 0"<<endl;
		//cout<<"#------------------------------------------"<<endl;
		if (BesselJ0_Zeros[0]/tLarge>uLim)
		{
			//cout<<"# in first if  lLim="<<lLim<<"  uLim="<<uLim<<endl;
			resultG=gaussianIntegrate_gsl(*this, lLim, uLim,100);
			//resultR=rombergIntegrate(*this,lLim, uLim,accuracyR);
		}
		else
		{
			arg=BesselJ0_Zeros[0]/tLarge;
			for(iL=0; lLim>(BesselJ0_Zeros[iL]/tLarge); iL++);
			for(iS=0; lLim>(BesselJ0_Zeros[iS]/tSmall);iS++);
			if((BesselJ0_Zeros[iS]/tSmall)< (BesselJ0_Zeros[iL]/tLarge))
				arg=BesselJ0_Zeros[iS]/tSmall;
			else
				arg=BesselJ0_Zeros[iL]/tLarge;
			//cout<<"#argument="<<arg<<"    iL="<<iL<<"   iS="<<iS<<"   lower limit="<<lLim<<endl;
			//number result1=rombergIntegrate(*this, lLim, arg,accuracyR,10E-8);
			//cout<<"# the first arg is="<<arg<<endl;
			//resultR+=rombergIntegrate(*this, lLim, arg,accuracyR);
			//cout<<"   resultG="<<resultG<<endl;

			//cout<<"#  lLim="<<lLim<<"   arg="<<arg<<endl;
			number result1=gaussianIntegrate_gsl(*this, lLim, arg,200);
			//int counter=0;
			resultG+=result1;
			//cout<<0<<"\t"<<resultG<<endl;
			while((arg<uLim)&&(iL<100000)&&(iS<100000)&&!terminate)
			{
				number lowerArg=arg;
				//cout<<"#----------------------------"<<endl;
				//cout<<"#iL="<<iL<<endl;
				if(tLarge==tSmall)
				{
					//counter++;
					arg=BesselJ0_Zeros[iL++]/tLarge;
					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<endl;
					number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<endl;
					//cout<<"#result1="<<result1<<endl;
					accuracyG=accuracyG-20;
					if(accuracyG<20)
						accuracyG=20;
					
					//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
					resultG+=result1;
					//cout<<counter<<"\t"<<resultG<<endl;
					if(abs(result1/resultG)<epsPP && lowerArg<arg)
					{
						//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
						terminate=true;
					}
					//resultR+=result2;
				}
				else
				{
					//counter++;
					if((BesselJ0_Zeros[iL]/tLarge)<(BesselJ0_Zeros[iS]/tSmall))
					{
						arg=BesselJ0_Zeros[iL]/tLarge;
						//cout<<"#     "<<BesselJ0_Zeros[iL]/tLarge<<"      "<<BesselJ0_Zeros[iS]/tSmall<<endl;
						//cout<<"#arg="<<arg<<endl;
						iL++;
					}
					else
					{
						arg=BesselJ0_Zeros[iS]/tSmall;
						//cout<<"#in else arg="<<arg<<endl;
						iS++;
					}

					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<endl;
					number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
					
					accuracyG=accuracyG-20;
					if(accuracyG<20)
						accuracyG=20;
					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<endl;
					//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<"   resultG="<<result1<<"   resultR="<<result2<<endl;
					resultG+=result1;
					//cout<<counter<<"\t"<<resultG<<endl;
					if(abs(result1/resultG)<epsPP && lowerArg<arg)
					{
						//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
						terminate=true;
					}
					//resultR+=result2;
				}
				//cout<<"# lowerArg="<<lowerArg<<" arg="<<arg<<endl;
			}

			if(arg>uLim)
			{
				//cout<<"# uLim="<<uLim<<"  arg="<<arg<<endl;
				resultG-=gaussianIntegrate_gsl(*this,uLim, arg,accuracyG);
				//cout<<counter+1<<"\t"<<resultG<<endl;
			}
			//cout<<"#argument="<<arg<<"  iL="<<iL<<"  iS="<<iS<<"  uLim="<<uLim<<endl;
// 			if (theta1==theta2)
// 			{
// 				//cout<<" # theta1=theta2, noiseTheta="<<noiseTheta<<endl;
// 				resultG+=noiseTheta;
// 			}
		}
	}

	else if((orderB1==0)&&(orderB2==4))
	{
		number epsPM=1e-7;
		if(theta1==theta2)
			epsPM=1e-12/(theta1/arcmin);
 		//cout<<"#epsPM="<<epsPM<<endl;
 		//cout<<"#theta1="<<theta1/arcmin<<"   theta2="<<theta2/arcmin<<endl;
 		//cout<<"#lLim="<<lLim<<"   uLim="<<uLim<<endl;
 		//cout<<"#--------------------------------------------------------------------"<<endl;
		int i0=0;
		int i4=0;
		//int counter=0;
		if(BesselJ0_Zeros[0]/theta1<BesselJ4_Zeros[0]/theta2)
			arg=BesselJ0_Zeros[0]/theta1;
		else
			arg=BesselJ4_Zeros[0]/theta2;
		//cout<<"# case 0"<<endl;
		if (arg>uLim)
		{
			cout<<"# in the first if"<<endl;
			resultG=gaussianIntegrate_gsl(*this, lLim, uLim,100);
			//resultR=rombergIntegrate(*this,lLim, uLim,accuracyR);
		}
		else
		{
			for(i0=0; lLim>(BesselJ0_Zeros[i0]/theta1); i0++);
			for(i4=0; lLim>(BesselJ4_Zeros[i4]/theta2); i4++);
			if((BesselJ4_Zeros[i4]/theta2)< (BesselJ0_Zeros[i0]/theta1))
				arg=BesselJ4_Zeros[i4]/theta2;
			else
				arg=BesselJ0_Zeros[i0]/theta1;
 			//cout<<"#lLim="<<lLim<<"    arg="<<arg<<endl;
			number result1=gaussianIntegrate_gsl(*this, lLim, arg,100);
			//number result2=rombergIntegrate(*this, lLim, arg,accuracyR);
			resultG+=result1;
 			//cout<<counter<<"\t"<<result1<<"\t"<<resultG<<endl;
			//resultR+=result2;
			//cout<<"#lLim="<<lLim<<"  arg="<<arg<<"  resultG="<<resultG<<endl;

			while((arg<uLim)&&(i0<100000)&&(i4<100000)&&!terminate)
			{
				//counter++;
				number lowerArg=arg;
				if((BesselJ0_Zeros[i0]/theta1)<(BesselJ4_Zeros[i4]/theta2))
				{
					arg=BesselJ0_Zeros[i0]/theta1;
					//cout<<"#     "<<BesselJ0_Zeros[iL]/tLarge<<"      "<<BesselJ0_Zeros[iS]/tSmall<<endl;
					//cout<<"#arg="<<arg<<endl;
					i0++;
				}
				else
				{
					arg=BesselJ4_Zeros[i4]/theta2;
					//cout<<"#in else arg="<<arg<<endl;
					i4++;
				}

 				//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<endl;
				number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
				accuracyG=accuracyG-20;
				if(accuracyG<20)
					accuracyG=20;
				if(abs(result1/resultG)<epsPM && lowerArg<arg)
				{
						//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
						terminate=true;
				}
				//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
				//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<"   resultG="<<result1<<"   resultR="<<result2<<endl;
				resultG+=result1;
 				//cout<<counter<<"\t"<<result1<<"\t"<<resultG<<endl;
				//resultR+=result2;
			}
		}
		if(arg>uLim)
		{
			//cout<<"#uLim="<<uLim<<"   arg="<<arg<<endl;
			resultG-=gaussianIntegrate_gsl(*this, uLim, arg,accuracyG);
		}
		//cout<<"#argument="<<arg<<"    i0="<<i0<<"   i4="<<i4<<"   upper limit="<<uLim<<endl;
	}

	else if ((orderB1==4)&&(orderB2==0))
	{
		number epsPM=1e-7;
		if(theta1==theta2)
			epsPM=1e-12/(theta1/arcmin);
		int i0=0;
		int i4=0;
		if(BesselJ0_Zeros[0]/theta2<BesselJ4_Zeros[0]/theta1)
			arg=BesselJ0_Zeros[0]/theta2;
		else
			arg=BesselJ4_Zeros[0]/theta1;
		//cout<<"# case 0"<<endl;
		if (arg>uLim)
		{
			resultG=gaussianIntegrate_gsl(*this, lLim, uLim,100);
		}
		else
		{
			for(i0=0; lLim>(BesselJ0_Zeros[i0]/theta2); i0++);
			for(i4=0; lLim>(BesselJ4_Zeros[i4]/theta1);i4++);
			if((BesselJ4_Zeros[i4]/theta1)< (BesselJ0_Zeros[i0]/theta2))
				arg=BesselJ4_Zeros[i4]/theta1;
			else
				arg=BesselJ0_Zeros[i0]/theta2;
			//cout<<"#argument="<<arg<<"    iL="<<iL<<"   iS="<<iS<<"   lower limit="<<lLim<<endl;
			resultG+=gaussianIntegrate_gsl(*this, lLim, arg,100);
			//resultR+=rombergIntegrate(*this, lLim, arg,accuracyR);
			//cout<<"#lower arg="<<lLim<<"   arg="<<arg<<endl;
			//cout<<"#resultR="<<resultR<<"   resultG="<<resultG<<endl;

			while((arg<uLim)&&(i0<100000)&&(i4<100000)&&!terminate)
			{
				number lowerArg=arg;
				if((BesselJ0_Zeros[i0]/theta2)<(BesselJ4_Zeros[i4]/theta1))
				{
					arg=BesselJ0_Zeros[i0]/theta2;
					i0++;
				}
				else
				{
					arg=BesselJ4_Zeros[i4]/theta1;
					i4++;
				}

				number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
				accuracyG=accuracyG-20;
				if(accuracyG<20)
					accuracyG=20;
				//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
				resultG+=result1;
				if(abs(result1/resultG)<epsPM && lowerArg<arg)
				{
					//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
					terminate=true;
				}
				//resultR+=result2;
				//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<"   resultG="<<result1<<"   resultR="<<result2<<"  difference="<<result1-result2<<endl;
				//resultR+=rombergIntegrate(*this, lowerArg, arg,accuracyR);
			}
		}
		if(arg>uLim)
			resultG-=gaussianIntegrate_gsl(*this, uLim, arg,accuracyG);
		//cout<<"#argument="<<arg<<"    i0="<<i0<<"   i4="<<i4<<"   upper limit="<<uLim<<endl;
		//cout<<theta1<<"\t"<<theta2<<"\t"<<(resultG-resultR)/resultG<<endl;
	}


	else if((orderB1==4)&&(orderB2==4))
	{
		number epsMM=1e-8;
// 		cout<<"#----------------------------------------------------------------------"<<endl;
// 		cout<<"#----------------------------------------------------------------------"<<endl;
// 		cout<<"#----------------------------------------------------------------------"<<endl;
// 		cout<<"#----------------------------------------------------------------------"<<endl;
// 		cout<<"#theta1="<<theta1/arcmin<<"  theta2="<<theta2/arcmin<<endl;
		accuracyG=20;
		int iL=0;
		int iS=0;
		number tLarge=0.;
		number tSmall=0.;
		if(theta1>theta2)
		{
			tLarge=theta1;
			tSmall=theta2;
		}
		else
		{
			tLarge=theta2;
			tSmall=theta1;
		}
		epsMM=epsMM*(tSmall/tLarge)*(tSmall/tLarge);
// 		cout<<"#epsMM="<<epsMM<<endl;
		//cout<<"#tLarge="<<tLarge<<"  tSmall="<<tSmall<<endl;
// 		cout<<"# case 4"<<endl;
		//int counter=0;
		if (BesselJ4_Zeros[0]/tLarge>uLim)
		{
			resultG=gaussianIntegrate_gsl(*this, lLim, uLim,100);
			//resultR=rombergIntegrate(*this,lLim, uLim,accuracyR);
		}
		else
		{
			arg=BesselJ4_Zeros[0]/tLarge;
			for(iL=0; lLim>(BesselJ4_Zeros[iL]/tLarge); iL++);
			for(iS=0; lLim>(BesselJ4_Zeros[iS]/tSmall);iS++);
			if((BesselJ4_Zeros[iS]/tSmall)< (BesselJ4_Zeros[iL]/tLarge))
				arg=BesselJ4_Zeros[iS]/tSmall;
			else
				arg=BesselJ4_Zeros[iL]/tLarge;
			//cout<<"#argument="<<arg<<"    iL="<<iL<<"   iS="<<iS<<"   lower limit="<<lLim<<endl;
			//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 				cout<<"#lLim="<<lLim<<"    arg="<<arg<<endl;
			

			resultG+=gaussianIntegrate_gsl(*this, lLim, arg,60);
			//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 				cout<<0<<"\t"<<resultG<<endl;
			//resultR+=rombergIntegrate(*this, lLim, arg,accuracyR);
			//cout<<"#resultR="<<resultR<<"   resultG="<<result<<endl;

			while((arg<uLim)&&(iL<100000)&&(iS<100000)&&!terminate)
			{
				//counter++;
				number lowerArg=arg;
				//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 					cout<<"#iL="<<iL<<"#iS="<<iS<<endl;
				if(tLarge==tSmall)
				{
					//cout<<"#equal theta"<<endl;
					arg=BesselJ4_Zeros[iL++]/tLarge;
					number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
					//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 						cout<<"#lowerArg="<<lowerArg<<"  arg="<<arg<<endl;
					//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
					resultG+=result1;
					//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 						cout<<counter<<"\t"<<resultG<<endl;
					if(abs(result1/resultG)<epsMM && lowerArg<arg )
					{
						//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
						terminate=true;
					}
					//resultR+=result2;
				}
				else
				{
					//counter++;
					if((BesselJ4_Zeros[iL]/tLarge)<(BesselJ4_Zeros[iS]/tSmall))
					{
						arg=BesselJ4_Zeros[iL]/tLarge;
						//cout<<"#     "<<BesselJ0_Zeros[iL]/tLarge<<"      "<<BesselJ0_Zeros[iS]/tSmall<<endl;
						//cout<<"#arg="<<arg<<endl;
						iL++;
					}
					else
					{
						arg=BesselJ4_Zeros[iS]/tSmall;
						//cout<<"#in else arg="<<arg<<endl;
						iS++;
					}

					//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 						cout<<"#lowerArg="<<lowerArg<<"  arg="<<arg<<endl;
					number result1=gaussianIntegrate_gsl(*this, lowerArg, arg,accuracyG);
					
// 					accuracyG=accuracyG-20;
// 					if(accuracyG<20)
// 						accuracyG=20;
					//number result2=rombergIntegrate(*this, lowerArg, arg,accuracyR);
					//cout<<"#lowerArg="<<lowerArg<<"   arg="<<arg<<"   resultG="<<result1<<"   resultR="<<result2<<endl;
					resultG+=result1;
					//if(theta1/arcmin==1. && theta2/arcmin==400.0)
// 						cout<<counter<<"\t"<<resultG<<endl;
					if(abs(result1/resultG)<epsMM && lowerArg<arg)
					{
						//cout<<"#result1="<<result1<<"  resultG="<<resultG<<endl;
						terminate=true;
					}
					//resultR+=result2;
				}
			}
			//cout<<"#argument="<<arg<<"    iL="<<iL<<"   iS="<<iS<<"   upper limit="<<uLim<<endl;
			if(arg>uLim)
			{
				resultG-=gaussianIntegrate_gsl(*this,uLim,arg,accuracyG);
// 				cout<<"#uLim="<<uLim<<"  arg="<<arg<<endl;
			}
// 			if (theta1==theta2)
// 			{
// 				//cout<<"theta1=theta2, noise="<<noise<<endl;
// 				resultG+=noiseTheta;
// 			}
			//resultR-=rombergIntegrate(*this,uLim,arg,accuracyG);
			//cout<<"#final resultG="<<resultG<<"  final resultR="<<resultR<<endl<<endl;
		}
	}

	else
	{
		clog<<"orderB1 or orderB2 not valid"<<endl;
		exit(1);
	}
	//writeInteg=false;
	//cout<<"# OrderB1="<<orderB1<<"   OrderB2="<<orderB2<<endl;
	//cout<<" # resultG="<<resultG<<endl;
	//cout<<"#A="<<A<<endl;
	//cout<<"#resultG/A/2./pi="<<resultG/A/2./pi<<endl;
	return resultG/A/2./pi;
}

// void KsiCovariance::calKsi()
// {
// 	clog<<"calculating Ksi in KsiCovariance"<<endl;
	
// 	KsiPlus Kp();
// 	KsiMinus Km();
// 	if(HistRedshift)
// 	{
// 		clog<<"HistRedshift"<<endl;
// 		Kp.setredshiftAndNoiseHist(pofz_mat_vec);
// 		Km.setredshiftAndNoiseHist(pofz_mat_vec);
// 	}
// 	Ksi_m_vec.clear();
// 	Ksi_p_vec.clear();
// 	int r=0;
// 	for(int bin1=0;bin1<nBins;bin1++)
// 	{
// 		for(int bin2=bin1;bin2<nBins;bin2++)
// 		{
// 			Ksi_p_vec.push_back(Kp);
// 			Ksi_m_vec.push_back(Km);
// 			Ksi_p_vec[r].MakeKsi(bin1,bin2);
//  			Ksi_m_vec[r].MakeKsi(bin1,bin2);
// 			r++;
// 		}
// 	}
// }


matrix KsiCovariance::CalCov(bool Logscale1)
{
	MixedTerm=false;
	matrix CMT(2*nKsiBins*nPairs,2*nKsiBins*nPairs);
	Logscale=Logscale1;
	clog<<"calculating the covariance in KsiCovariance"<<endl;
	if(checkPower())
	{
		if(Logscale)
		{
			clog<<"log binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);

							number theta1=thetamin;

							for(int i=2*nKsiBins*p1,n=0;i<nKsiBins*(2*p1+1);i++,n++)
							{
								// if(logbinning)
								// {
								// 	number theta_binned=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+0.5));
								// 	// 							clog<<"logbining, theta_binned="<<theta_binned/arcmin<<endl;
								// 	theta_max=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+1.));
								// 	theta_min=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta));
								// }
								// else
								// {
								// 	number theta_binned=(tmax-tmin)/(nTheta_binned)*(itheta+0.5)+tmin;
								// 	// 							clog<<"linbinning, theta_binned="<<theta_binned/arcmin<<endl;
								// 	theta_max=theta_binned+(tmax-tmin)/(nTheta_binned)*0.5;
								// 	theta_min=theta_binned-(tmax-tmin)/(nTheta_binned)*0.5;
								// }
								
								theta1=exp(log(thetamin)+LogCoef*(n+0.5));
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=2*nKsiBins*p2+n,m=n; j<nKsiBins*(2*p2+1); j++,m++)
								{
									theta2=exp(log(thetamin)+LogCoef*(m+0.5));
									cout<<"#theta2="<<theta2/arcmin<<endl;
									//string CovFileName=
									//		folderName+string("/Covs/Cov_")
									//		+toString(theta1/arcmin,8)+string("-")
									//		+toString(theta2/arcmin,8)+string(".bin");
									matrix Cov(2,2);
									//ifstream Covfile(CovFileName.c_str());
									//if(Covfile.fail())
									//{
										Cov.load(0,0,value(theta1,theta2,0,0));
										Cov.load(1,1,value(theta1,theta2,4,4));
										Cov.load(0,1,value(theta1,theta2,0,4));
										Cov.load(1,0,value(theta1,theta2,4,0));
									//	Cov.store(CovFileName.c_str());
									//}
									//else
									//	Cov.restore(CovFileName.c_str());
									//Covfile.close();

									CMT.load(j,i,Cov.get(0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));

		///very important: the way the covariance is built is very important if +- is replaced by -+ 
		//it is no longer positive defenite
									
									CMT.load(j+nKsiBins,i,Cov.get(0,1));
									CMT.load(j,i+nKsiBins,Cov.get(1,0));
									CMT.load(j+nKsiBins-m+n,i-n+m,CMT.get(j,i+nKsiBins));
									CMT.load(i,j+nKsiBins,CMT.get(j+nKsiBins,i));

									CMT.load(j+nKsiBins,i+nKsiBins,Cov.get(1,1));
									CMT.load(i+nKsiBins,j+nKsiBins,CMT.get(j+nKsiBins,i+nKsiBins));
									CMT.load(i+nKsiBins+m-n,j+nKsiBins+n-m,CMT.get(i+nKsiBins,j+nKsiBins));
									CMT.load(j+nKsiBins-m+n,i+nKsiBins-n+m,CMT.get(j+nKsiBins,i+nKsiBins));
								}
							}
						}
					}
				}
			}
		}
		else
		{
			clog<<"Linear binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							number theta1=thetamin;

							for(int i=2*nKsiBins*p1,n=0;i<nKsiBins*(2*p1+1);i++,n++)
							{
								theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=2*nKsiBins*p2+n,m=n; j<nKsiBins*(2*p2+1); j++,m++)
								{
									theta2=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*m;

									cout<<"#theta1="<<theta1/arcmin<<"   theta2="<<theta2/arcmin<<endl;
									// string CovFileName=
									// 		folderName+string("/Covs2/Cov_")
									// 		+toString(theta1/arcmin,8)+string("-")
									// 		+toString(theta2/arcmin,8)+string(".bin");
									matrix Cov(2,2);
									// ifstream Covfile(CovFileName.c_str());
									// if(Covfile.fail())
									// {
										Cov.load(0,0,value(theta1,theta2,0,0));
										Cov.load(1,1,value(theta1,theta2,4,4));
										Cov.load(0,1,value(theta1,theta2,0,4));
										Cov.load(1,0,value(theta1,theta2,4,0));
										// Cov.store(CovFileName.c_str());
									// }
									// else
										// Cov.restore(CovFileName.c_str());
									// Covfile.close();

									CMT.load(j,i,Cov.get(0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));

		///very important: the way the covariance is built is very important if +- is replaced by -+ it is no longer positive defenite
									CMT.load(j+nKsiBins,i,Cov.get(0,1));
									CMT.load(j,i+nKsiBins,Cov.get(1,0));
									CMT.load(j+nKsiBins-m+n,i-n+m,CMT.get(j,i+nKsiBins));
									CMT.load(i,j+nKsiBins,CMT.get(j+nKsiBins,i));

									CMT.load(j+nKsiBins,i+nKsiBins,Cov.get(1,1));
									CMT.load(i+nKsiBins,j+nKsiBins,CMT.get(j+nKsiBins,i+nKsiBins));
									CMT.load(i+nKsiBins+m-n,j+nKsiBins+n-m,CMT.get(i+nKsiBins,j+nKsiBins));
									CMT.load(j+nKsiBins-m+n,i+nKsiBins-n+m,CMT.get(j+nKsiBins,i+nKsiBins));
								}
							}
						}
					}
				}
			}
		}

		if(!NoNoise)
		{
			if(Logscale)
			{
				for(int bin1=0; bin1<nBins; bin1++)
				{
					for(int bin2=bin1; bin2<nBins; bin2++)
					{
						for(int bin3=bin1; bin3<nBins; bin3++)
						{
							for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
							{
								
								int n,m;
								deltaNoise=delta(bin1,bin3)*delta(bin2,bin4)
											+delta(bin1,bin4)*delta(bin2,bin3);
								///This is not right, need to change use COSEBIs covXiB
								///The shot noise term should be:
								//sigma_e^4/(2 pi theta DeltaTheta A nbar^2) x deltanoise/2
								number noiseTheta1=deltaNoise*noise_vec[bin1]*noise_vec[bin2]/A/pi;
								//the pair considered for the Eparam_vec
								int p1=calP(nBins,bin1,bin2);
								int p2=calP(nBins,bin3,bin4);
								number theta1=thetamin;
								for(int i=2*nKsiBins*p1,n=0;i<nKsiBins*(2*p1+1);i++,n++)
								{
									theta1=exp(log(thetamin)+LogCoef*(n+0.5));
									number deltaTheta=exp(log(thetamin)+LogCoef*(n+1.))-exp(log(thetamin)+LogCoef*n);
									number noiseTheta=noiseTheta1/deltaTheta/theta1;
									//cout<<"#theta1="<<theta1/arcmin<<endl;
									number Cov=CMT.get(i,i)+noiseTheta;
									CMT.load(i,i,Cov);
									Cov=CMT.get(i+nKsiBins,i+nKsiBins)+noiseTheta;
									CMT.load(i+nKsiBins,i+nKsiBins,Cov);
								}
							}
						}
					}
				}
			}
			else
			{
				number deltaTheta=(thetamax-thetamin)/(nKsiBins-1.0);
				for(int bin1=0; bin1<nBins; bin1++)
				{
					for(int bin2=bin1; bin2<nBins; bin2++)
					{
						for(int bin3=bin1; bin3<nBins; bin3++)
						{
							for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
							{
								
								int n,m;
								deltaNoise=delta(bin1,bin3)*delta(bin2,bin4)
											+delta(bin1,bin4)*delta(bin2,bin3);
								number noiseTheta1=deltaNoise*noise*noise/A/pi;
								//the pair considered for the Eparam_vec
								int p1=calP(nBins,bin1,bin2);
								int p2=calP(nBins,bin3,bin4);
								number theta1=thetamin;
								for(int i=2*nKsiBins*p1,n=0;i<nKsiBins*(2*p1+1);i++,n++)
								{
									theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
									number noiseTheta=noiseTheta1/theta1/deltaTheta;
									cout<<"#noiseTheta="<<noiseTheta<<endl;
									cout<<"#theta1="<<theta1/arcmin<<endl;
									number Cov=CMT.get(i,i)+noiseTheta;
									CMT.load(i,i,Cov);
									Cov=CMT.get(i+nKsiBins,i+nKsiBins)+noiseTheta;
									CMT.load(i+nKsiBins,i+nKsiBins,Cov);
								}
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



matrix KsiCovariance::CalCovPlus(bool Logscale1)
{
	MixedTerm=false;
	matrix CMT(nKsiBins*nPairs,nKsiBins*nPairs);
	Logscale=Logscale1;
	clog<<"calculating the covariance in KsiCovariance"<<endl;

	if(checkPower())
	{
		if(Logscale)
		{
			clog<<"log binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);

							number theta1=thetamin;

							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
							{
								// if(logbinning)
								// {
								// 	number theta_binned=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+0.5));
								// 	// 							clog<<"logbining, theta_binned="<<theta_binned/arcmin<<endl;
								// 	theta_max=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta+1.));
								// 	theta_min=exp(log(tmin)+log(tmax/tmin)/(nTheta_binned)*(itheta));
								// }
								// else
								// {
								// 	number theta_binned=(tmax-tmin)/(nTheta_binned)*(itheta+0.5)+tmin;
								// 	// 							clog<<"linbinning, theta_binned="<<theta_binned/arcmin<<endl;
								// 	theta_max=theta_binned+(tmax-tmin)/(nTheta_binned)*0.5;
								// 	theta_min=theta_binned-(tmax-tmin)/(nTheta_binned)*0.5;
								// }
								
								theta1=exp(log(thetamin)+LogCoef*(n+0.5));
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
								{
									theta2=exp(log(thetamin)+LogCoef*(m+0.5));
									cout<<"#theta2="<<theta2/arcmin<<endl;
									//string CovFileName=
									//		folderName+string("/Covs/Cov_")
									//		+toString(theta1/arcmin,8)+string("-")
									//		+toString(theta2/arcmin,8)+string(".bin");
									//matrix Cov(2,2);
									//ifstream Covfile(CovFileName.c_str());
									//if(Covfile.fail())
									//{
										//Cov.load(0,0,value(theta1,theta2,0,0));
										// Cov.load(1,1,value(theta1,theta2,4,4));
										// Cov.load(0,1,value(theta1,theta2,0,4));
										// Cov.load(1,0,value(theta1,theta2,4,0));
									//	Cov.store(CovFileName.c_str());
									//}
									//else
									//	Cov.restore(CovFileName.c_str());
									//Covfile.close();

									CMT.load(j,i,value(theta1,theta2,0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));

		///very important: the way the covariance is built is very important if +- is replaced by -+ 
		//it is no longer positive defenite
									
									// CMT.load(j+nKsiBins,i,Cov.get(0,1));
									// CMT.load(j,i+nKsiBins,Cov.get(1,0));
									// CMT.load(j+nKsiBins-m+n,i-n+m,CMT.get(j,i+nKsiBins));
									// CMT.load(i,j+nKsiBins,CMT.get(j+nKsiBins,i));

									// CMT.load(j+nKsiBins,i+nKsiBins,Cov.get(1,1));
									// CMT.load(i+nKsiBins,j+nKsiBins,CMT.get(j+nKsiBins,i+nKsiBins));
									// CMT.load(i+nKsiBins+m-n,j+nKsiBins+n-m,CMT.get(i+nKsiBins,j+nKsiBins));
									// CMT.load(j+nKsiBins-m+n,i+nKsiBins-n+m,CMT.get(j+nKsiBins,i+nKsiBins));
								}
							}
						}
					}
				}
			}
		}
		else
		{
			clog<<"Linear binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							number theta1=thetamin;

							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
							{
								theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
								{
									theta2=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*m;

									cout<<"#theta1="<<theta1/arcmin<<"   theta2="<<theta2/arcmin<<endl;
									// string CovFileName=
									// 		folderName+string("/Covs2/Cov_")
									// 		+toString(theta1/arcmin,8)+string("-")
									// 		+toString(theta2/arcmin,8)+string(".bin");
									//matrix Cov(2,2);
									// ifstream Covfile(CovFileName.c_str());
									// if(Covfile.fail())
									// {
										//Cov.load(0,0,value(theta1,theta2,0,0));
										//Cov.load(1,1,value(theta1,theta2,4,4));
										//Cov.load(0,1,value(theta1,theta2,0,4));
										//Cov.load(1,0,value(theta1,theta2,4,0));
										// Cov.store(CovFileName.c_str());
									// }
									// else
										// Cov.restore(CovFileName.c_str());
									// Covfile.close();

									CMT.load(j,i,value(theta1,theta2,0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));

		///very important: the way the covariance is built is very important if +- is replaced by -+ it is no longer positive defenite
									// CMT.load(j+nKsiBins,i,Cov.get(0,1));
									// CMT.load(j,i+nKsiBins,Cov.get(1,0));
									// CMT.load(j+nKsiBins-m+n,i-n+m,CMT.get(j,i+nKsiBins));
									// CMT.load(i,j+nKsiBins,CMT.get(j+nKsiBins,i));

									// CMT.load(j+nKsiBins,i+nKsiBins,Cov.get(1,1));
									// CMT.load(i+nKsiBins,j+nKsiBins,CMT.get(j+nKsiBins,i+nKsiBins));
									// CMT.load(i+nKsiBins+m-n,j+nKsiBins+n-m,CMT.get(i+nKsiBins,j+nKsiBins));
									// CMT.load(j+nKsiBins-m+n,i+nKsiBins-n+m,CMT.get(j+nKsiBins,i+nKsiBins));
								}
							}
						}
					}
				}
			}
		}


		if(!NoNoise)
		{
			if(Logscale)
			{
				for(int bin1=0; bin1<nBins; bin1++)
				{
					for(int bin2=bin1; bin2<nBins; bin2++)
					{
						for(int bin3=bin1; bin3<nBins; bin3++)
						{
							for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
							{
								
								int n,m;
								deltaNoise=delta(bin1,bin3)*delta(bin2,bin4)
											+delta(bin1,bin4)*delta(bin2,bin3);
								//This is not right, need to change use COSEBIs covXiB
								number noiseTheta1=deltaNoise*noise_vec[bin1]/A/pi;
								//the pair considered for the Eparam_vec
								int p1=calP(nBins,bin1,bin2);
								int p2=calP(nBins,bin3,bin4);
								number theta1=thetamin;
								for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
								{
									theta1=exp(log(thetamin)+LogCoef*(n+0.5));
									number deltaTheta=exp(log(thetamin)+LogCoef*(n+1.))-exp(log(thetamin)+LogCoef*n);
									number noiseTheta=noiseTheta1/deltaTheta/theta1;
									cout<<"#theta1="<<theta1/arcmin<<endl;
									number Cov=CMT.get(i,i)+noiseTheta;
									CMT.load(i,i,Cov);
									Cov=CMT.get(i+nKsiBins,i+nKsiBins)+noiseTheta;
									CMT.load(i+nKsiBins,i+nKsiBins,Cov);
								}
							}
						}
					}
				}
			}
			else
			{
				number deltaTheta=(thetamax-thetamin)/(nKsiBins-1.0);
				for(int bin1=0; bin1<nBins; bin1++)
				{
					for(int bin2=bin1; bin2<nBins; bin2++)
					{
						for(int bin3=bin1; bin3<nBins; bin3++)
						{
							for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
							{
								
								int n,m;
								deltaNoise=delta(bin1,bin3)*delta(bin2,bin4)
											+delta(bin1,bin4)*delta(bin2,bin3);
								number noiseTheta1=deltaNoise*noise*noise/A/pi;
								//the pair considered for the Eparam_vec
								int p1=calP(nBins,bin1,bin2);
								int p2=calP(nBins,bin3,bin4);
								number theta1=thetamin;
								for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
								{
									theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
									number noiseTheta=noiseTheta1/theta1/deltaTheta;
									cout<<"#noiseTheta="<<noiseTheta<<endl;
									cout<<"#theta1="<<theta1/arcmin<<endl;
									number Cov=CMT.get(i,i)+noiseTheta;
									CMT.load(i,i,Cov);
									Cov=CMT.get(i+nKsiBins,i+nKsiBins)+noiseTheta;
									CMT.load(i+nKsiBins,i+nKsiBins,Cov);
								}
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

matrix KsiCovariance::CalCovPlusMixed(bool Logscale1)
{
	MixedTerm=true;
	matrix CMT(nKsiBins*nPairs,nKsiBins*nPairs);
	Logscale=Logscale1;
	clog<<"calculating the mixed term of covariance of Xip in KsiCovariance"<<endl;

	if(checkPower())
	{
		if(Logscale)
		{
			clog<<"log binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);

							number theta1=thetamin;

							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
							{
								theta1=exp(log(thetamin)+LogCoef*(n+0.5));
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
								{
									theta2=exp(log(thetamin)+LogCoef*(m+0.5));
									cout<<"#theta2="<<theta2/arcmin<<endl;
									CMT.load(j,i,value(theta1,theta2,0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));
								}
							}
						}
					}
				}
			}
		}
		else
		{
			clog<<"Linear binning for theta"<<endl;
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							int n,m;
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							number theta1=thetamin;

							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
							{
								theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
								cout<<"#theta1="<<theta1/arcmin<<endl;
								number theta2=thetamin;
								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
								{
									theta2=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*m;

									cout<<"#theta1="<<theta1/arcmin<<"   theta2="<<theta2/arcmin<<endl;
									CMT.load(j,i,value(theta1,theta2,0,0));
									CMT.load(i,j,CMT.get(j,i));
									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
									CMT.load(j-m+n,i-n+m,CMT.get(j,i));
								}
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


bool KsiCovariance::checkPower()
{
	if(powerspectrum_vec.size())
		return true;
	return false;
}


void KsiCovariance::readNpairs(vector<string> FileName_vec,int nColumns)
{
	clog<<"in readNpairs in KsiCovariance"<<endl;
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



vector<vector<number> > KsiCovariance::readKsi(string FileName,int nColumns)
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
       			 //clog<<temp<<endl;
        		cols_temp++;
      		}
      		if(cols_temp>0)
      		{
      			//clog<<"cols_temp="<<cols_temp<<endl;
				Ksi_vecvec.push_back(Ksi_vec);
				nRow++;
			}
		}
	}
	//clog<<"nRows="<<nRow<<endl;
	//exit(1);
	KsiFile.close();
	return Ksi_vecvec;
}


// vector<vector<number> > KsiCovariance::readKsi(string FileName,int nColumns)
// {
// 	int nRow=0;
// 	string str;
// 	number temp;
// 	vector<vector<number> > Ksi_vecvec;
// 	ifstream KsiFile((FileName).c_str());
//   	clog<<"reading, "<<FileName<<endl;
// 	if(KsiFile.fail())
// 	{
// 		clog<<"error occured during opening: "<<FileName<<endl;
// 		exit(1);
// 	}

// 	string line;
// 	while(getline(KsiFile,line))
// 	{
// 		stringstream stream(line);
// 		str=stream.peek();
// 		if (str=="#")
// 		{
// 			//nColumns++;
// // 			string line;
// // 			getline(KsiFile,line);
// // 			clog<<line<<endl;
// 		}
// 		else
// 		{
// 			int cols_temp=0;
//       		number temp;
//       		vector<number> Ksi_vec;
//       		while(stream>>temp)
//       		{
//         		Ksi_vec.push_back(temp);
//        			 //clog<<temp<<endl;
//         		cols_temp++;
//       		}
//       		if(cols_temp>0)
//       		{
//       			//clog<<"cols_temp="<<cols_temp<<endl;
// 				Ksi_vecvec.push_back(Ksi_vec);
// 				nRow++;
// 			}
// 		}
// 	}
// 	//clog<<"nRows="<<nRow<<endl;
// 	//exit(1);
// 	KsiFile.close();
// 	return Ksi_vecvec;
// }

// vector<vector<number> > KsiCovariance::readfile(string FileName,int nColumns)
// {
// 	int nRow=0;
// 	string str;
// 	number temp;
// 	vector<vector<number> > w2ww_vecvec;
// 	ifstream w2wwFile((FileName).c_str());
// //  	clog<<"reading, "<<FileName<<endl;
// 	if(w2wwFile.fail())
// 	{
// 		clog<<"error occured during opening: "<<FileName<<endl;
// 		exit(1);
// 	}
// 	///lets change this
// 	string line;
// 	while(getline(w2wwFile,line))
// 	{
// 		stringstream stream(line);
// 		str=stream.peek();
// 		if (str=="#")
// 		{
// 			//nColumns++;
// // 			string line;
// // 			getline(KsiFile,line);
// // 			clog<<line<<endl;
// 		}
// 		else
// 		{
// 			vector<number> w2ww_vec(nColumns);
// 			for(int nCol=0; nCol<nColumns; nCol++)
// 			{
// 				stream>>temp;
// //  				cout<<temp<<"\t";
// 				w2ww_vec[nCol]=temp;
// 			}
// // 			cout<<endl;
// // 			exit(1);
// 			w2ww_vecvec.push_back(w2ww_vec);
// 			nRow++;
// 		}
// 	}
// // 	clog<<"read Ksifile"<<endl;
// 	w2wwFile.close();
// // 	Ksi_vecvec.pop_back();
// //  	clog<<"Ksi_vecvec.size()="<<Ksi_vecvec.size()<<"  Ksi_vecvec[0].size()="<<Ksi_vecvec[0].size()<<endl;
// 	return w2ww_vecvec;
// }

///This doesn't work for w2ww yet there are two sets of thetas: theta12 and theta13 and both need to be in range
vector<int> KsiCovariance::FindMinMaxTheta(vector<vector<number> > Ksi_vec)
{
// 	int row=0;
// 	clog<<"finding min and max theta that are closest to thetamin and thetamax"<<endl;
	vector<int> index(2);
	clog.precision(10);
	number h=Ksi_vec[1][0]-Ksi_vec[0][0];
 	clog<<"h="<<h<<endl;
	number tmin=thetamin;
	number tmax=thetamax;
	if(Ksi_vec[0][0]>tmin || Ksi_vec[Ksi_vec.size()-1][0]<tmax)
	{
		clog<<"the theta range is not sufficient:"
			<<Ksi_vec[0][0]<<"   "<<Ksi_vec[Ksi_vec.size()-1][0]<<endl;
		exit(1);
	}
	//clog<<"Ksi_vec[900][0]="<<Ksi_vec[900][0]<<endl;
	//clog<<"Ksi_vec[901][0]="<<Ksi_vec[901][0]<<endl;
	
	clog<<"tmin="<<tmin<<endl;
	clog<<"tmax="<<tmax<<endl;
	int ntmin=ceil((tmin-Ksi_vec[0][0])/h);
	int ntmax=Ksi_vec.size()-1;
	clog<<"at the start ntmin="<<ntmin<<endl;
	clog<<"tmax-Ksi_vec[0][0]="<<tmax-Ksi_vec[0][0]<<endl;
	//clog<<"at the start ntmax="<<ntmax<<endl;
	while(tmin>Ksi_vec[ntmin][0])
		ntmin++;
	clog<<"Ksi_vec["<<ntmin<<"]="<<Ksi_vec[ntmin][0]<<endl;
	while(tmax<Ksi_vec[ntmax][0])
		ntmax--;
	//ntmax--;
	//clog<<"ntmax="<<ntmax<<endl;
	clog<<"Ksi_vec["<<ntmax<<"]="<<Ksi_vec[ntmax][0]<<endl;
	index[0]=ntmin;
	index[1]=ntmax;
	exit(1);
	return index;
}



void KsiCovariance::setIntegLimitTheta()
{
	int ilow=0;
	int ihigh=Theta_vec.size()-1;
	if(Theta_vec[ihigh]<thetamax)
	{
		clog<<"max Theta_vec is smaller than thetamax, exiting now"<<endl;
		exit(1);
	}
	if(Theta_vec[ilow]>thetamin)
	{
		clog<<"min Theta_vec is larger than thetamin, exiting now"<<endl;
		exit(1);
	}
	while(Theta_vec[ilow]<thetamin)
		ilow++;
	while(Theta_vec[ihigh]>thetamax)
		ihigh++;
	itheta_low=ilow;
	N_b_range=ihigh-ilow+1;
	integ_limit_theta.clear();
	integ_limit_theta.push_back(itheta_low);
	int itheta=itheta_low;
	for(int iksibin=1;iksibin<=nKsiBins; iksibin++)
	{
		number theta_edge=0.;
		if(Logscale)
			theta_edge=exp(log(thetamin)+LogCoef*iksibin);
		else
			theta_edge=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*iksibin;
		while(Theta_vec[itheta]<theta_edge)
			itheta++;
		integ_limit_theta.push_back(itheta);
	}
}

// ///value of covariance
// number KsiCovariance::valueInputw2ww(int itheta1,int itheta2,char WhichKsi)
// {
// 	number result=0.;
// // make info for triangle bins
//   //   number dlna = log(flmax/flmin)/bins;
//   //   number dlnb = log(flmax/flmin)/bins;
//   //   number dc   = 2.*pi/bins;
    
//   //   for(int n1=0;n1<bins;n1++)
// 		// for(int n2=0;n2<bins;n2++)
// 		//     for(int n3=0;n3<bins;n3++)
// 		//     {
// 		// 		number a = .5*(exp(log(flmin)+dlna*(n1+1.))+exp(log(flmin)+dlna*n1));
// 		// 		number b = .5*(exp(log(flmin)+dlnb*(n2+1.))+exp(log(flmin)+dlnb*n2));
// 		// 		number c  = dc*(n3+.5);
				
// 		// 		long index = (long) (n1+n2*bins+n3*bins*bins);
				
// 		// 		number da = exp(log(flmin)+dlna*(n1+1.))-exp(log(flmin)+dlna*n1);
// 		// 		number db = exp(log(flmin)+dlnb*(n2+1.))-exp(log(flmin)+dlnb*n2);
				
// 		// 		W123.load(0,index,a); theta12
// 		// 		W123.load(1,index,b); theta13
// 		// 		W123.load(2,index,c); phi3
// 		// 		W123.load(3,index,da);
// 		// 		W123.load(4,index,db);
// 		// 		W123.load(5,index,dc);
// 		//     }
// 	///algorithm for this: 
// 	// We have a vertor of all thetas that are in the w2ww file: Theta_vec
// 	// We also have a integ_limit_theta that knows theta indices in Theta_vec are in the range for
// 	// the set of theta_lo and theta_high which are the same for theta12 and theta13
// 	// and we know the Number of bins, N_b, in the w2ww file.
// 	// and we have rearranged w2ww to be like this for N_b=2:
// 	// theta1 theta1 phi1
// 	// theta1 theta1 phi2
// 	// theta1 theta2 phi1
// 	// theta1 theta2 phi2
// 	// theta2 theta1 phi1
// 	// theta2 theta1 phi2
// 	// theta2 theta2 phi1
// 	// theta2 theta2 phi2
// 	//
// 	// Then we go through the w2ww matrix till we find a theta that is larger than theta1_low
// 	// Then we go through the seocnd column of w2ww till we find a theta that is higher than theta2_low
// 	// There exactly N_b rows with this condition and they have variying phi values. 
// 	// Then we sum over all these N_b rows x xi(theta23)

//  	number integ_theta12=0.;
//  	for(int itheta_12=integ_limit_theta[itheta1]; itheta_12<integ_limit_theta[itheta1+1]; itheta_12++)
// 	{
// 		number theta12=Theta_vec[itheta_12];
// 		///this is the starting row 
// 		int row12=itheta_12*N_b*N_b;
// 		clog<<"theta12="<<w2ww_mat_vec[redshiftPair].get(0,row12)<<endl;
// 	 	number integ_theta13=0.;
// 		for(int itheta_13=integ_limit_theta[itheta2]; itheta_13<integ_limit_theta[itheta2+1]; itheta_13++)
// 		{
// 			number theta13=Theta_vec[itheta_13];
// 			///this is the starting row 
// 			int row1213=row12+ (itheta_13*N_b);
// 			clog<<"theta13="<<w2ww_mat_vec[redshiftPair].get(1,row1213)<<endl;
// 			clog<<"theta13-1="<<w2ww_mat_vec[redshiftPair].get(1,row1213-1)<<endl;
// 			number integ_sumphi=0.;
// 			for(int ibin=0; ibin<N_b; ibin++)
// 			{
// 				number phi3=w2ww_mat_vec[redshiftPair].get(2,row1213+ibin);
// 				number theta23=sqrt(theta12*theta12+theta13*theta13-2.*theta12*theta13*cos(phi3));
// 				number ksi=0.;
// 				if(WhichKsi=='p')
// 				{
// 					ksi=Ksi_p_vec[redshiftPair].value(theta23);
// 				}
// 				else if(WhichKsi=='m')
// 				{
// 					ksi=Ksi_m_vec[redshiftPair].value(theta23);
// 				}
// 				number integ= w2ww_mat_vec[redshiftPair].get(6,row1213+ibin)*ksi;
// 				integ_sumphi+=integ;
// 				clog<<"theta12="<<w2ww_mat_vec[redshiftPair].get(0,row1213+ibin)
// 					<<" theta13="<<w2ww_mat_vec[redshiftPair].get(1,row1213+ibin)
// 					<<" phi3="<<w2ww_mat_vec[redshiftPair].get(2,row1213+ibin)<<endl;
// 			}
// 			integ_theta13+=integ_sumphi;
// 		}
// 		integ_theta12+=integ_theta13;
// 	}
// 	result=integ_theta12/(Np[itheta1]*Np[itheta2]);
// 	return result;
// }

// void KsiCovariance::readw2ww(vector<string> FileName_vec,int nColumns)
// {
// 	clog<<"in readw2ww in ksiCovariance"<<endl;
// 	clog<<"FileName_vec.size is:"<<FileName_vec.size()<<endl;
// 	nPairs=FileName_vec.size();
// 	w2ww_mat_vec.clear();
// 	for(int r=0; r<nPairs; r++)
// 	{
// 		vector<vector<number> > w2ww=readfile(FileName_vec[r],nColumns);
// 		matrix w2ww_mat(7,w2ww.size()-1);
// 		for(int i=0;i<w2ww.size(); i++)
// 		{
// 			w2ww_mat.load(0,i,w2ww[i][0]);
// 			w2ww_mat.load(1,i,w2ww[i][1]);
// 			w2ww_mat.load(2,i,w2ww[i][2]);
// 			w2ww_mat.load(3,i,w2ww[i][3]);
// 			w2ww_mat.load(4,i,w2ww[i][4]);
// 			w2ww_mat.load(5,i,w2ww[i][5]);
// 			w2ww_mat.load(5,i,w2ww[i][6]);
// 			w2ww_mat.load(5,i,w2ww[i][7]);
// 		}
// 		w2ww_mat_vec.push_back(w2ww_mat);
// 		clog<<"w2ww_mat[0]="<<w2ww_mat.get(0,0)<<"   w2ww_mat[1]="<<w2ww_mat.get(1,0)<<endl;
// 	}
// 	N_b=int(pow(w2ww_mat_vec[0].rows,1./3.));
// 	for(int i=0;i<N_b; i++)
// 	{
// 		Theta_vec[i]=w2ww_mat_vec[0].get(1,i*N_b);
// 		clog<<"Theta_vec["<<i<<"]="<<Theta_vec[i];
// 	}
// }

// ///This is not done yet
// matrix KsiCovariance::CalCovPlusMixed_Inputw2ww(bool Logscale1)
// {
// 	MixedTerm=true;
// 	setKsi(sigma8,omega_matter,omega_lambda,w0,nindex,h100,omega_baryon);


// 	matrix CMT(nKsiBins*nPairs,nKsiBins*nPairs);
// 	Logscale=Logscale1;
// 	clog<<"calculating the mixed term of covariance of Xip in KsiCovariance"<<endl;
// 	setIntegLimitTheta();

// 	string Binning;
// 	if (Logscale)
// 		Binning="LogBins-";
// 	else
// 		Binning="LinearBins-";

// 	// string CMTFileName=
// 	// 		folderName+string("/CMTPlusMixed")+Binning
// 	// 		+powerswitchName+string("-")+survey+("-")+toString(nKsiBins)+string("_")
// 	// 		+toString(begin,2)+string("-")+toString(end,2)+string("_")
// 	// 		+toString(nBins)+string("-")
// 	// 		+toString(sigma8,4)+string("-")+toString(omega_matter,4)+string("-")
// 	// 		+toString(omega_lambda,4)+string("-")+toString(w0,4)+string("-")
// 	// 		+toString(nindex,4)+string("-")+toString(h100,4)+string("-")
// 	// 		+toString(omega_baryon,4)+string("-")+toString(sigmaE,4)
// 	// 		+string("_")+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2)
// 	// 		+string(".bin");


// 	// clog<<"in calcov, file name is: "<<CMTFileName<<endl;
// 	// ifstream CMTfile(CMTFileName.c_str());
// 	// if(CMTfile.fail())
// 	// {
// 		// clog<<" CMT file failed"<<endl;
// 		// calPower();
// 		if(Logscale)
// 		{
// 			clog<<"log binning for theta"<<endl;
// 			for(int bin1=0; bin1<nBins; bin1++)
// 			{
// 				for(int bin2=bin1; bin2<nBins; bin2++)
// 				{
// 					for(int bin3=bin1; bin3<nBins; bin3++)
// 					{
// 						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
// 						{
// 							int n,m;
// 							rp1=calP(nBins,bin1,bin3);
// 							rp2=calP(nBins,bin2,bin4);
// 							rp3=calP(nBins,bin1,bin4);
// 							rp4=calP(nBins,bin2,bin3);
// 							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
// 							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
// 							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
// 							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
// 							int p1=calP(nBins,bin1,bin2);
// 							int p2=calP(nBins,bin3,bin4);

// 							number theta1=thetamin;

// 							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
// 							{
// 								int itheta1=n;
// 								// theta1=exp(log(thetamin)+LogCoef*(n+0.5));
// 								cout<<"#theta1="<<theta1/arcmin<<endl;
// 								number theta2=thetamin;
// 								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
// 								{
// 									int itheta2=m;
// 									theta2=exp(log(thetamin)+LogCoef*(m+0.5));
// 									cout<<"#theta2="<<theta2/arcmin<<endl;
// 									CMT.load(j,i,valueInputw2ww(itheta1,itheta2,'p'));
// 									CMT.load(i,j,CMT.get(j,i));
// 									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
// 									CMT.load(j-m+n,i-n+m,CMT.get(j,i));
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		else
// 		{
// 			clog<<"Linear binning for theta"<<endl;
// 			for(int bin1=0; bin1<nBins; bin1++)
// 			{
// 				for(int bin2=bin1; bin2<nBins; bin2++)
// 				{
// 					for(int bin3=bin1; bin3<nBins; bin3++)
// 					{
// 						for(int  bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
// 						{
// 							int n,m;
// 							rp1=calP(nBins,bin1,bin3);
// 							rp2=calP(nBins,bin2,bin4);
// 							rp3=calP(nBins,bin1,bin4);
// 							rp4=calP(nBins,bin2,bin3);
// 							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
// 							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
// 							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
// 							delta4noise=delta(bin2,bin3)*noise_vec[bin2];
// 							int p1=calP(nBins,bin1,bin2);
// 							int p2=calP(nBins,bin3,bin4);
// 							number theta1=thetamin;

// 							for(int i=nKsiBins*p1,n=0;i<nKsiBins*(p1+1);i++,n++)
// 							{
// 								theta1=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*n;
// 								cout<<"#theta1="<<theta1/arcmin<<endl;
// 								number theta2=thetamin;
// 								int itheta1=n;
// 								for(int j=nKsiBins*p2+n,m=n; j<nKsiBins*(p2+1); j++,m++)
// 								{
// 									theta2=thetamin+(thetamax-thetamin)/(nKsiBins-1.)*m;
// 									int itheta2=m;
// 									cout<<"#theta1="<<theta1/arcmin<<"   theta2="<<theta2/arcmin<<endl;
// 									CMT.load(j,i,valueInputw2ww(itheta1,itheta2,'p'));
// 									CMT.load(i,j,CMT.get(j,i));
// 									CMT.load(i+m-n,j+n-m,CMT.get(i,j));
// 									CMT.load(j-m+n,i-n+m,CMT.get(j,i));
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	// 	CMT.store(CMTFileName.c_str());
// 	// }
// 	// else 
// 	// {
// 	// 	CMT.restore(CMTFileName.c_str());
// 	// }
// 	// CMTfile.close();
// 	return CMT;
// }







// matrix KsiCovariance::CovInterpolate(int order,int NewNKsiBins,bool Logscale)
// {
// 	///this is the input matrix which is going to be interpolated to find the desired Covariance
// 	///using a polynomail fit to the covariance
// 	matrix CMT=CalCov(Logscale);
// 	matrix CMTNEW(2*NewNKsiBins*nPairs,2*NewNKsiBins*nPairs);
// 	int Nbins=nKsiBins*nPairs;
// 	matrix CMTPlus(Nbins,Nbins);
// 	matrix CMTMinus(Nbins,Nbins);
// 	matrix CMTPM(Nbins,Nbins);
// 	//change this for tomography
// 	for(int i=0; i<Nbins; i++)
// 	{
// 		for(int j=0; j<Nbins; j++)
// 		{
// 			CMTPlus.load(i,j,CMT.get(i,j));
// 			CMTMinus.load(i,j,CMT.get(i+Nbins,j+Nbins));
// 			CMTPM.load(i,j,CMT.get(i,j+Nbins));
// 		}
// 	}
// 	CMTPlus.printOut((folderName+string("/CovPlus_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	CMTPlus.polyfit(order);
// 	CMTPlus.printOut((folderName+string("/CovPlusPoly_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	CMTMinus.printOut((folderName+string("/CovMinus_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	CMTMinus.polyfit(order);
// 	CMTMinus.printOut((folderName+string("/CovMinusPoly_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	CMTPM.printOut((folderName+string("/CovPM_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	CMTPM.polyfit(order);
// 	CMTPM.printOut((folderName+string("/CovPMPoly_")+toString(order)
// 			+string("_")+survey+string("-")
// 			+toString(nKsiBins)+string("_")
// 			+toString(nBins)+string("_")+toString(thetamin/arcmin,2)
// 			+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
// 	return CMTNEW;
// }


int KsiCovariance::calP(int nBins,int fbin,int sbin)
{
	int p=fbin*nBins;
	for(int i=0; i<fbin; i++)
	{
		p-=i;
	}
	return p+sbin-fbin;
}
