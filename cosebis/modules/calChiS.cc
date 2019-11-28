#include "calChiS.h"

number calChiS(matrix theory,matrix data,matrix Cov)
{
	if(data.rows!=theory.rows || data.columns!=theory.columns)
	{
		clog<<"!!!!!!!!!!!!!!!!!!!data and theory size don't match!!!!!!!!!!!!!!!!!"<<endl;
	}
	// else
	// {
	// 	clog<<"!!!!!!!!!!!!!!!!!!They match!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	// }
	matrix Delta=data-theory;
	matrix iCov=Cov.inverse();
	matrix chiS=Delta.t()*iCov*Delta;
	return chiS.get(0);
}
