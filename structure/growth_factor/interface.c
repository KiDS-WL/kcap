#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include "growthfactor.h"


//Module to calculate the linear growth factor D, and linear growth rate, f. Where D, f are defined by the growth of a
//linear perturbation, delta, with scale factor a: 
//delta(a') = delta(a)*(D(a')/D(a)) and f = dlnD/dlna
//Anyone using Komatsu's CRL library should note: growth_factor_crl = D *(1+z) and growth_rate_crl = f/(1+z)


const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * like = LIKELIHOODS_SECTION;
const char * growthparameters = GROWTH_PARAMETERS_SECTION;

typedef struct growth_config {
	double zmin;
	double zmax;
	double dz;
	int nz;
} growth_config;

void reverse(double * x, int n)
{
	for (int i=0; i<n/2; i++){
		double tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i] = tmp;
	}
}

growth_config * setup(c_datablock * options)
{
	int status = 0;
	growth_config * config = malloc(sizeof(growth_config));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmin", 0.0, &(config->zmin));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax", 3.0, &(config->zmax));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "dz", 0.01, &(config->dz));
	config->nz = (int)((config->zmax-config->zmin)/config->dz)+1;
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);

	printf("Will calculate f(z) and d(z) in %d bins (%lf:%lf:%lf)\n", config->nz, config->zmin, config->zmax, config->dz);
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
        if (status){
                fprintf(stderr, "Please specify the redshift in the growth function module.\n");
                exit(status);
        }
	return config;

} 
    
int execute(c_datablock * block, growth_config * config)
{

	int i,status=0;
	double w,wa,omega_m,omega_v;
	int nz = config->nz;
	

	//read cosmological params from datablock
    status |= c_datablock_get_double_default(block, cosmo, "w", -1.0, &w);
    status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0, &wa);
    status |= c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
    status |= c_datablock_get_double_default(block, cosmo, "omega_lambda", 1-omega_m, &omega_v);

	if (status){
		fprintf(stderr, "Could not get required parameters for growth function (%d)\n", status);
		return status;
	}

	//allocate memory for single D, f and arrays as function of z
	double *a = malloc(nz*sizeof(double));
	double *dz = malloc(nz*sizeof(double));
	double *fz = malloc(nz*sizeof(double));
	double *z = malloc(nz*sizeof(double));


	// output D and f over a range of z
	// we do this in reverse becase the growth code is expecting
	// increasing values of "a"
	// we reverse afterwards
	for (i=0;i<nz;i++){
		z[nz-1-i] = config->zmin + i*config->dz;
		a[nz-1-i] = 1.0/(1+z[nz-1-i]);
	}


	status = get_growthfactor(nz, a, omega_m, omega_v, w, wa, dz, fz);
	
	reverse(a,nz);
	reverse(z,nz);
	reverse(dz,nz);
	reverse(fz,nz);


	status |= c_datablock_put_double_array_1d(block,growthparameters, "d_z", dz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "f_z", fz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "z", z, nz);
	free(fz);
	free(dz);
	free(z);
	free(a);

return status;
}


int cleanup(growth_config * config)
{
	free(config);
	return 0;
}
