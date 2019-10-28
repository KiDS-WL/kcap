#include "nrutil.h"

//----- the following is borrowed from the numerical recipes

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	cerr << "Numerical Recipes run-time error..."<<endl;
        cerr << error_text << endl;
	cerr << "...now exiting to system...\n"<<endl;
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

number *Vector(long nl, long nh)
/* allocate a number vector with subscript range v[nl..nh] */
{
	number *v;

	v=(number *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(number)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(number *v, long nl, long nh)
/* free a number vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

number **Matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	number **m;

	/* allocate pointers to rows */
	m=(number **) malloc((size_t)((nrow+NR_END)*sizeof(number*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(number *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(number)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_Matrix(number **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
