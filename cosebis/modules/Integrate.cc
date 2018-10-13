#include "Integrate.h"

//gsl Gaussian legendre integration
number gaussianIntegrate_gsl(function_cosebis& f,number a, number b, int N)
{
  if (a==b) return 0.0;

  number xi;
  number wi;
  number result = 0.0;
  gsl_integration_glfixed_table * t= gsl_integration_glfixed_table_alloc(N);
  for(int m=0;m<N;m++)
  {
    gsl_integration_glfixed_point(a, b, m, &xi, &wi, t);
    result+=wi*f.integrant(xi);
    //cout<<m+1<<'\t'<<xi<<'\t'<<wi<<endl;
  }
  gsl_integration_glfixed_table_free(t);
  return result;
}