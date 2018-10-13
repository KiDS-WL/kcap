#include "unistd.h"
#include "matrix.h"

 //--------- trigonometric functions

 number _cos(const number a) {return cos(a);}

 number _sin(const number a) {return sin(a);}

 number _tan(const number a) {return tan(a);}

 number _acos(const number a){return acos(a);}

 number _asin(const number a){return asin(a);}

 number _atan(const number a){return atan(a);}

 //----------------------------------------------------------------------
void matrix_error(char error_text[])
{
  cerr << "Error"<<endl;
  cerr << error_text << endl;
  cerr << ", Closing now\n"<<endl;
  exit(1);
}

 int alphaToUnsigned(number x)
 {
  return (int) (x*RESOLUTION/(2.0*pi));
 }

int betaToUnsigned(number x)
{
  return (int) ((x+pi/2.0)*RESOLUTION/pi);
}
 
int rotate(int alpha,  int beta,
                 int euler1, int euler2, int euler3, 
                 bool deinstall)
{
  static int *lookuptable=NULL;
  static bool wasHere = false;

  if (deinstall)
  {
   //--- use deinstall=true to free allocated memory !

   if (lookuptable==NULL) delete[] lookuptable;

   wasHere = false;

   cerr << "deinitialize lookuptable: rotate" << endl;

   return 0;
  }

  if (wasHere)
  {
   int alpha1, beta1, ab2;
   int alpha2, beta2, ab3;
   int alpha3, beta3;

   //--- use buffer to calculate rotation fast

   //-- convert to int according to resolution
  
   alpha1 = (alpha+euler1)%RESOLUTION;
   beta1  = beta;

   ab2    = lookuptable[alpha1+beta1*RESOLUTION];
   alpha2 = ab2&(RESOLUTION-1);
   beta2  = ab2-alpha2;
   alpha2 = (alpha2+euler2)%RESOLUTION;

   ab3    = lookuptable[alpha2+beta2*RESOLUTION];
   alpha3 = ab3&(RESOLUTION-1);
   beta3  = ab3-alpha3;
   alpha3 = (alpha3+euler3)%RESOLUTION;

   //--- returns direction (alpha,beta) as compact int
   //--- = alpha+beta*RESOLUTION
   //--- take care to convert your texture to this format and you
   //--- can use this as texture index right away !

   return beta3+alpha3;
  }

  clog << "initializing lookuptable: rotate" << endl;
  
  //--- initialize buffer on first call

  lookuptable = new int[RESOLUTION*RESOLUTION];

  int alpha_int;
  int beta_int;
  int alpha_int2;
  int beta_int2;
  number    alpha1;
  number    alpha2;
  number    beta1;
  number    beta2;
  number    x,y,z;

  for(alpha_int=0;alpha_int<RESOLUTION;alpha_int++)
  for(beta_int =0;beta_int <RESOLUTION;beta_int++)
  {
   //--- direction in spherical rads

   alpha1 = alpha_int*(2.0*pi)/RESOLUTION;
   beta1  = beta_int*pi/RESOLUTION-pi/2.0;

   //--- direction in cartesian coords (rotated about x-axis)

   x      = cos(alpha1)*cos(beta1);
   y      = sin(beta1);
   z      = sin(alpha1)*cos(beta1);

   //--- spherical rads in rotated coordinate system
    
   beta2  = asin(z);
   alpha2 = acos(x/cos(beta2));
   if (y/cos(beta2)<0) alpha2 = 2.0*pi-alpha2;

   //--- and converted to integer according to resolution

   alpha_int2 = alphaToUnsigned(alpha2);
   beta_int2  = betaToUnsigned(beta2);

   //--- write to "transfer"-matrix
 
   lookuptable[alpha_int+beta_int*RESOLUTION] = 
    alpha_int2+beta_int2*RESOLUTION;     
  }   
  
  wasHere = true;
  
  return 0;
}

 /*
 
   ----------------------------------------------------------------------
   ------ class matrix --------------------------------------------------
   ----------------------------------------------------------------------
  
*/

ostream& operator<<(ostream& s, matrix& a)
{
  a.printOut(s);
 
  return s;
} 

matrix operator+(const matrix& a, number n)
{
  matrix c=a;
 
  c.add(n);

  return c;
}

matrix operator+(number n, const matrix& a)
{
  matrix c=a;
 
  c.add(n);

  return c;
 }

matrix operator-(const matrix& a, number n)
{
  matrix c=a;
 
  c.sub(n);

  return c;
}

matrix operator-(number n, const matrix& a)
{
  matrix c=a;
 
  c.sub(n);

  return c;
}

matrix operator^(const matrix& a,const matrix& b)
{
  matrix c=rot3d(a);

  c.mul(b);

  return c;
 }

matrix operator^(const matrix& a,int n)
{
  matrix   c=a;
  int i;

  for(i=1;i<n;i++) c.mul(a);

  return c;
 }  

 matrix operator+(const matrix& a, const matrix& b)
 { 
  matrix c=a;

  c.add(b);

  return c;
 }

 matrix operator-(const matrix& a, const matrix& b)
 {  
  matrix c=a;

  c.sub(b);

  return c;
 }

 matrix operator*(const matrix& a, const matrix& b)
 {
  matrix c=a;

  //--- standard product for matrices

  if (a.columns==b.rows)
  {
   c.mul(b);

   return c;
  }

  //--- yields standard scalar product of two vectors
  //--- and an generalization for matrices of same number of columns/rows

  if (a.columns==b.columns || a.rows==b.rows)
  {
   c.transpose();
   c.mul(b);

   return c;
  }

  //--- returns simply a on failure

  c.error("incompatible matrices: mul");

  return a;
 }

 matrix operator*(const matrix& a, number n)
 {
  matrix c=a;

  c.mul(n);

  return c;
 }

 matrix operator*(number n, const matrix& a)
 {
  matrix c=a;

  c.mul(n);

  return c;
 }


 matrix operator/(const matrix& a, const matrix& b)
 {
  matrix c=a;
  matrix d=b;

  c.mul(d.inverse());

  return c;
 }

 matrix operator/(number n, const matrix& a)
 {
  matrix c=a;
  
  c.inverse();
  c.mul(n);

  return c;
 }

 matrix operator/(const matrix& a, number n)
 {
  matrix c=a;

  c.div(n);

  return c;
 }

 bool operator==(matrix& a, const matrix& b)
 {
  return a.equals(b);
 }

 //-------- matrix methods ----------------------------------------------

 number* matrix::getElements() {return elements;}

 number matrix::scalar()
 {
  if (columns==1 && rows==1) return elements[0];
   else return mod();
 }

 void matrix::error(const char *message)
 {
  cerr << message << endl;

  exit(1);
 }

 matrix::matrix(int x, int y)
 {
  columns = x;
  rows    = y;

  elements= new number[x*y];
  neutral(); 
 }

matrix::matrix(const char* file)
{
  columns = 0;
  rows    = 0;
  elements= NULL;
  restore(file);
}

void matrix::resize(int x,int y)
{
  clear();

  columns = x;
  rows    = y;

  elements= new number[x*y];

  neutral();
}

void matrix::resize(int x)
{
  resize(1,x);
  neutral();
}

matrix::matrix(int y)
 {
  columns = 1;
  rows    = y;

  elements = new number[y];
  neutral();
 }

matrix::matrix()
 {
  columns = 0;
  rows    = 0;
  
  elements = NULL;
 }

matrix::matrix(const matrix& a)
 {
  columns  = a.columns;
  rows     = a.rows;
  elements = new number[rows*columns];
  
  memcpy(elements,a.elements,sizeof(number)*columns*rows);
 }

 void matrix::operator=(const matrix& a)
 {
  if (elements!=NULL) delete[] elements;

  columns = a.columns;
  rows    = a.rows;
  elements= new number[columns*rows];

  memcpy(elements,a.elements,sizeof(number)*columns*rows);
 }

matrix::~matrix() 
 {
   clear(); 
 }

void matrix::clear()
 {
  if (elements!=NULL) delete[] elements;

  elements = NULL;
  columns  = 0;
  rows     = 0;
 }

 void matrix::load(number *newElements)
 {
   memcpy(elements,newElements,sizeof(number)*columns*rows);
 }

 void matrix::load(int x,int y,number value)
 {
    // if (x>=columns)
    //   error("too many columns in matrix");
    // if (y>=rows)
    //   error("too many rows in matrix");
    elements[inArray(x,y)] = value;
 }

 void matrix::load(int x,number value)
 {
   // //clog<<"value"<<endl;
   // if (x>=columns)
   //    error("too many columns in matrix");
   load(x,0,value);
 }  

 void matrix::operator=(number *newElements)
 {
    load(newElements);
 }

int matrix::inArray(int x,int y) const
{
  if (!x && x>=columns || !y && y>=rows)
  {
    cerr << "matrix: invalid reference to element: "
	    << " x="  << x
	    << " y=" << y
	    << endl;

    x%=columns;
    y%=rows;
  }

  return x+y*columns;
}

number matrix::get(int x) const
{
    return get(x,0);
}

number matrix::get(int x,int y) const
{
  return elements[inArray(x,y)];
}

void matrix::add(const matrix& b)
{
  int n;

  if (b.columns*b.rows==columns*rows)
    for(n=0;n<columns*rows;n++) 
      elements[n]+= b.elements[n];
  else error("incompatible matrices: add");

}

void matrix::add(number n)
{
  matrix c(columns,rows);

  add(c*n);
}

 void matrix::add(int x,int y,number n)
 {
   elements[x+y*columns]+=n; 
 }

 void matrix::add(int x,number n)
 {
   elements[x]+=n;
 }

 void matrix::sub(const matrix& b)
 {
  int n;

 if (b.columns*b.rows==columns*rows)
   for(n=0;n<columns*rows;n++) 
     elements[n]-= b.elements[n];
 else error("incompatible matrices: sub");  

 }

 void matrix::sub(number n)
 {
  matrix c(columns,rows);

  sub(c*n);
 }

 void matrix::sub(int x,int y,number n)
 {
   elements[x+y*columns]-=n;
 }

 void matrix::sub(int x,number n)
 {
   elements[x]-=n;
 }

 void matrix::mul(const matrix& b)
 {
  if (b.rows!=columns)
   error("incompatible matrices: mul");
  
  number*  newElements = new number[b.columns*rows];
  long     col,row,i;
  number   result;

  for(row=0;row<(long)rows;row++)
    for(col=0;col<(long)b.columns;col++)
    {
     for(result=0.,i=0;i<(long)columns;i++)
       result+= 
         elements[i+columns*row]*b.elements[col+b.columns*i];   

     newElements[col+row*b.columns] = result;
    }

  columns = b.columns;
  if (elements!=NULL) 
    delete[] elements;

  elements = newElements;
 }

 void matrix::mul(number f)
 {
  int n;

  for(n=0;n<columns*rows;n++) 
    elements[n]*= f;
 }

 void matrix::mul(int x,int y,number n)
 {
   elements[x+y*columns]*=n;
 }

 void matrix::mul(int x,number n)
 {
   elements[x]*=n;
 }

 void matrix::div(number f)
 {
  int n;

  for(n=0;n<columns*rows;n++) 
    elements[n]/= f;
 }

 void matrix::div(int x,int y,number n)
 {
   elements[x+y*columns]/=n;
 }

 void matrix::div(int x,number n)
 {
   elements[x]/=n;
 }

 bool matrix::equals(const matrix& b)
 {
  if (columns!=b.columns || rows!=b.rows) return false;

  int n;
  for(n=0;n<columns*rows;n++) 
    if (elements[n]!=b.elements[n]) return false;

  return true;
 }

 void matrix::transpose()
 {
  number *newElements = new number[rows*columns];

  int n;
  int m;

  for(n=0;n<rows;n++)
    for(m=0;m<columns;m++)  
      newElements[n+rows*m] = 
	elements[m+n*columns];

  memcpy(elements,newElements,sizeof(number)*columns*rows);
  delete[] newElements;

  n       = rows;
  rows    = columns;
  columns = n;
 }

 void matrix::zero(number fill)
 {
  int n;

  for(n=0;n<rows*columns;n++) 
    elements[n] = fill;
 }

 void matrix::neutral()
 {
   long m,n;
   for(n=0;n<(long)rows;n++)
     for(m=0;m<(long)columns;m++)
       elements[m+n*columns] = (m==n ? 1. : 0.);
 }  

 matrix matrix::symmetric()
 {
  matrix c(columns,rows);
  
  if (columns!=rows) 
    return c;

  long m,n;
  for(n=0;n<(long)rows;n++)
    for(m=0;m<(long)columns;m++)
      c.elements[m+n*columns] = 
	elements[m+n*columns]+elements[n+m*columns];
  
  c.mul(_half);

  return c;
 }

 matrix matrix::asymmetric()
 {
  matrix c(columns,rows);
  
  if (columns!=rows) return c;

  long m,n;
  for(n=0;n<(long)rows;n++)
    for(m=0;m<(long)columns;m++)
      c.elements[m+n*columns] = 
	elements[m+n*columns]-elements[n+m*columns];

  c.mul(_half);

  return c;
 }
  
 number matrix::mod()
 {
  number result = _zero;

  int n; 
  for(n=0;n<columns*rows;n++) 
    result = _add(result,_sqr(elements[n]));

  return _sqrt(result);
 }
/*
 matrix matrix::norm()
 {
  matrix c(columns,rows);
  c.load(elements);
  c.normalize();

  return c;
 }*/
/*
 void matrix::normalize()
 {
  number modulus = mod();

  if (modulus) 
    div(modulus);
 }*/

number  matrix::trace()
{
  number trace=0.;

  if (rows!=columns)
    return 0.;

  for(int n=0;n<rows;n++)
    trace += get(n,n);

  return trace;
}

void matrix::strip(int d)
{
  int n,m;
  
  for(n=0;n<rows;n++)
    for(m=0;m<columns;m++)
      if (abs(m-n)>d)
	load(m,n,0.);
}

 matrix matrix::t()
 {
  matrix b(columns,rows);
  b.load(elements);

  b.transpose();

  return b;
 }

 number matrix::max()
 {
   number maximum=0.;
   number value=0.;
   int n;
   int m;

   for(n=0;n<rows;n++)
     for(m=0;m<columns;m++)
       if (!n && !m)
	 maximum = get(m,n);
       else if ((value=get(m,n))>maximum) maximum = value;

   return maximum;
 }


 number matrix::min()
 {
   number minimum=0.;
   number value=0.;
   int n;
   int m;

   for(n=0;n<rows;n++)
     for(m=0;m<columns;m++)
       if (!n && !m)
	 minimum = get(m,n);
       else if ((value=get(m,n))<minimum) minimum = value;

   return minimum;
 }

 void  matrix::operator+=(const matrix& a)
 {
  matrix b = *this+a;

  *this=b;
 }

 void  matrix::operator-=(const matrix& a)
 {
  matrix b = *this-a;

  *this=b;
 }
   
 void  matrix::operator+=(number a)
 {
  matrix b = *this+a;
 
  *this=b;
 }

 void  matrix::operator*=(const matrix& a)
 {
  matrix b = *this*a;

  *this = b;
 }

 void  matrix::operator*=(number a)
 {
  matrix b = *this*a;

  *this = b;
 }
 /*
 void  matrix::operator/=(const matrix& a)
 {
  matrix b = *this/a;

  *this = b;
 }

 void  matrix::operator/=(number a)
 {
  matrix b = *this/a;

  *this = b;
 }
*/
void matrix::printOut(const char* filename,int prec)
{
  ofstream fhandler;
  fhandler.open(filename);
  if (!fhandler.fail())
    {
      printOut(fhandler,prec);
      fhandler.close();
    }
  else
    cerr << "matrix: could not open " << filename << " for writing" << endl;
}

void matrix::printOutNoNewLine(const char* filename,int prec)
{
  ofstream fhandler;
  fhandler.open(filename);
  if (!fhandler.fail())
    {
      printOutNoNewLine(fhandler,prec);
      fhandler.close();
    }
  else
    cerr << "matrix: could not open " << filename << " for writing" << endl;
}


void matrix::printOut(ostream& out, int prec)
{
  int n;
  int m;

  out.setf(ios::scientific | ios::showpos);
  out.precision(prec);
  
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    out << elements[m+n*columns] << "\t";
  out << endl;
}  
/*
//testing
void matrix::printOut(ostream& out, int prec)
{
  int n;
  int m;

  clog.setf(std::ios::scientific | std::ios::showpos);
  clog.precision(prec);
  clog<<"in printOut"<<endl;
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    clog << elements[m+n*columns] << "\t";
  clog << endl;
}  */

void matrix::printOutNoNewLine(ostream& out, int prec)
{
  int n;
  int m;

  out.setf(ios::scientific | ios::showpos | ios::floatfield);
  out.precision(prec);
  
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    out << elements[m+n*columns] << "\t";
  //out << endl;
}  


void matrix::store(const char* filename)
{
  FILE* fhandler = fopen(filename,"wb");

  if (fhandler==NULL)
  {
      cerr << "matrix: could not open filename " 
	         << filename<< endl;
      return;
  }

  fwrite((void*) &columns,sizeof(int),1,fhandler);
  fwrite((void*) &rows,sizeof(int),1,fhandler);
  
  number tmp;
  for(int n=0;n<columns*rows;n++)
  {
      tmp = (number) elements[n];
      fwrite((void*) &tmp,sizeof(number),1,fhandler);
  }

  fclose(fhandler);
}

void matrix::restore(const char* filename)
{
  clear();

  FILE* fhandler = fopen(filename,"rb");

  if (fhandler==NULL)
  { 
     clog<<"matrix: not able to restore: "<<filename<<endl;
	  exit(1);
     return;
  }

  fread((void*) &columns,sizeof(int),1,fhandler);
  fread((void*) &rows,sizeof(int),1,fhandler);

  elements = new number[columns*rows];
  number tmp;
  for(int n=0;n<columns*rows;n++)
  {
      fread((void*) &tmp,sizeof(number),1,fhandler);
      elements[n] = (number) tmp;
  }
      
  fclose(fhandler);
}

matrix matrix::inverse()
{
  matrix b(columns,rows);
  b.load(elements);

  b.invert();

  return b;
}

/*void matrix::storeFITS(string str)
{
  unlink(str.c_str());

  // to use FITS subroutines we have to store everything as number
  number* tmp = new number[columns*rows];
  for(int x=0;x<columns;x++)
    for(int y=0;y<rows;y++)
      tmp[x+y*columns] = (number)get(x,y);

  // make fits file
  fitsfile *fptr;
  long fpixel[2]= {1,1};
  long naxes[2] = {columns,rows};
  int  naxis=2, status=0;

  fits_create_file(&fptr,str.c_str(),&status);
  fits_create_img(fptr,FLOAT_IMG,naxis,naxes,&status);
  fits_write_pix(fptr,TFLOAT,fpixel,size(),(void*) tmp,&status);
  fits_close_file(fptr,&status);
  fits_report_error(stderr,status);
  delete[] tmp;
}*/

/*int matrix::restoreFITS(string str, int layer)
{
  clear();

  // copy & paste from a cfitsio website...
  fitsfile *fptr;
  long anaxes[3] = {1,1,1};
  int  status    = 0;
  int  anaxis;

  //Open the input file and create output file 
  fits_open_file(&fptr,str.c_str(), READONLY, &status);  
  fits_get_img_dim(fptr, &anaxis, &status);  //read dimensions 
  fits_get_img_size(fptr, 3, anaxes, &status);

  if (status)
    return 0;

  int npixels = anaxes[0]*anaxes[1]; 
  double* ptr = new double[npixels];
  if (anaxis>2)
  {
      // three-dimensional FITS file
      alarm(layer>=anaxes[2],"selected layer for matrix does not exit!");
      
      long fpixel[] = {1,1,layer+1};
      fits_read_pix(fptr, TDOUBLE, fpixel, npixels, NULL, ptr, NULL, &status);
  }
  else
  {
    // two-dimensional FITS file
    long fpixel[] = {1,1};
    int  hdutype;
    if (layer>1)
	       fits_movabs_hdu(fptr, layer, &hdutype, &status);

    if (status)
	     return 0;

    fits_read_pix(fptr, TDOUBLE, fpixel, npixels, NULL, ptr, NULL, &status);
  }
    
  long count = 0;
  resize(anaxes[0],anaxes[1]);
  for(int rows=0;rows<anaxes[1];rows++)
    for(int cols=0;cols<anaxes[0];cols++)
      load(cols,rows,ptr[count++]);
  fits_close_file(fptr, &status);

  delete[] ptr;

  // return the number of layers
  return anaxes[2];
}*/

/*
matrix matrix::inverse2()
{
  matrix b(columns,rows);
  b.load(elements);

  b.invert2();

  return b;
}

matrix matrix::inverse3(number thresh)
{
  matrix U,V,W;
  makesvd(*this,U,W,V);

   matrix Wmat(W.rows,W.rows);
  Wmat.zero(0.);
  for(int n=0;n<W.rows;n++)
    Wmat.load(n,n,W.get(n)>thresh?1./W.get(n):0.);

  return V*Wmat*U.t();
}

void matrix::getSVD(matrix& U,matrix& W,matrix& V)
{
  matrix Wvec;
  makesvd(*this,U,Wvec,V);
  W.resize(Wvec.rows,Wvec.rows);
  W.zero(0.);
  for(int n=0;n<Wvec.rows;n++)
    W.load(n,n,Wvec.get(n));  
}
*/
/*
void   matrix::getLU(matrix& L,matrix& U)
{
  L.resize(rows,rows);
  U.resize(rows,rows);
  L.zero();
  U.zero();

  matrix dummy = *this;
  dummy.LUdecompose();

  int m,n;
  for(n=0;n<rows;n++)
    for(m=0;m<rows;m++)
    {
    	if (m<=n)
    	  L.load(m,n,(m==n ? 1. : dummy.get(m,n)));
    	
    	if (m>=n)
    	  U.load(m,n,dummy.get(m,n));
    }
}
///This one has NR
void matrix::cholesky()
{
  if (!size())
    return;
  
  if (columns!=rows)
    error("Cholesky: matrix is not square matrix");

  int N = rows;

  number** a  = Matrix(1,N,1,N);
  number*  p  = Vector(1,N);
  
  int m,n;
  for(n=1;n<=N;n++)
    for(m=1;m<=N;m++)
      a[n][m] = get(m-1,n-1);

  choldc(a,N,p);
  zero();

  for(n=1;n<=N;n++)
    for(m=n;m<=N;m++)
      if (n!=m)
	load(m-1,n-1,a[m][n]);
      else
	load(m-1,n-1,p[n]);

  free_Matrix(a,1,N,1,N);
  free_vector(p,1,N);
}
//NR
number matrix::LUdecompose()
{
  if (!size())
    return 0.;
  
  if (columns!=rows)
    error("LUdecompose: matrix is not square matrix");
  
  int N = rows;

  number** a  = Matrix(1,N,1,N);
  int*    indx= new int[N+1];
  number  d;
  
  int m,n;
  for(n=1;n<=N;n++)
    for(m=1;m<=N;m++)
      a[n][m] = get(m-1,n-1);

  ludcmp(a,N,indx,&d);

  for(n=1;n<=N;n++)
    for(m=1;m<=N;m++)
      load(m-1,n-1,a[n][m]);

  free_Matrix(a,1,N,1,N);
  delete[] indx;

  return d;
}

matrix matrix::choleskydecomposed()
{
  matrix result = *this;

  result.cholesky();

  return result;
}

matrix matrix::LUdecomposed()
{
  matrix result = *this;

  result.LUdecompose();

  return result;
}

number  matrix::det()
{
  matrix dummy  = *this;
  number result = dummy.LUdecompose();

  int n;
  for(n=0;n<rows;n++)
      result*=dummy.get(n,n);
 
  return result;
}

matrix  matrix::diagonal()
{
  if (columns!=rows) 
    error("cannot diagonalize non-square matrix");

  matrix eigenvectors;
  matrix eigenvalues;

  eigenvalvect(*this,eigenvectors,eigenvalues);

  int m,n;
  for(n=0;n<rows;n++)
    for(m=0;m<columns;m++)
      load(m,n,(m==n ? eigenvalues.get(m) : 0.));

  return eigenvectors;
}

void matrix::invert3(number thresh)
{
  *this = this->inverse3(thresh);
}

void matrix::invert2()
{
  if (columns!=rows) 
    error("cannot invert non-square matrix");

  if (columns==1) 
  {
   elements[0] = _inverse(elements[0]);
   return;
  }
  
  int N = rows;
  number** a  = Matrix(1,N,1,N);
  int*    indx= new int[N+1];
  number* col = new number[N+1];
  number  d;
  
  int m,n;
  for(n=1;n<=N;n++)
    for(m=1;m<=N;m++)
      a[n][m] = get(m-1,n-1);

  ludcmp(a,N,indx,&d);

  int i,j;
  for(j=1;j<=N;j++)
  {
    for(i=1;i<=N;i++) 
	     col[i]=0.;
    col[j]=1.;
    lubksb(a,N,indx,col);
    for(i=1;i<=N;i++)
	     load(i-1,j-1,col[i]);
  }

  free_Matrix(a,1,N,1,N);
  delete[] indx;
  delete[] col;
 }
*/
 void matrix::invert()
 {
	if (columns!=rows) 
	error("cannot invert non-square matrix");

	if (columns==1) 
	{
		elements[0] = _inverse(elements[0]);
		return;
	}

	//--- invoke numerical recipes here

	number** matrix1 = Matrix(1,rows,1,rows);
	number** matrix2 = Matrix(1,rows,1,rows);

	long m,n;
	for(m=0;m<(long)rows;m++)
	for(n=0;n<(long)rows;n++)
	{
		matrix2[m+1][n+1] = (m==n ? 1.:0.);
		matrix1[m+1][n+1] = elements[n+m*rows];
	}

	gaussj(matrix1,rows,matrix2,rows);

	for(m=0;m<(long)rows;m++)
	for(n=0;n<(long)rows;n++)
	elements[n+m*rows] = matrix1[m+1][n+1];

	free_Matrix(matrix1,1,rows,1,rows);
	free_Matrix(matrix2,1,rows,1,rows);
 }

int matrix::size() const
{
  return columns*rows;
}

 matrix matrix::getRow(int row)
 {
   matrix result(columns,1);
   int n;

   for(n=0;n<columns;n++)
     result.load(n,row,get(n,row));

   return result;
 }
/*
number matrix::cond()
{
  if (!size())
    return 1.;

  matrix U,W,V;
  makesvd(*this,U,W,V);

  number min,max;
  for(int n=0;n<W.rows;n++)
    if (!n)
      min = max = W.get(n);
    else
      {
	min = (W.get(n)<min? W.get(n):min);
	max = (W.get(n)>max? W.get(n):max);
      }

  return max/min;
}*/

 matrix matrix::getColumn(int column)
 {
   matrix result(1,rows);
   int n;

   if (column>=columns)
     {
       cerr << "matrix: out of area column grab [STOPPED]!"
	    << endl;

       cin >> n;

       return result;
     }

   for(n=0;n<rows;n++)
     result.load(0,n,get(column,n));

   return result;
 }

void matrix::addColumnFromFile(const char* filename,int col,int cols)
{
  ifstream input;
  input.open(filename);
  if (input.fail())
    {
      cerr << "matrix: could not open file "
	   << filename
	   << endl;

      exit(1);
    }

  // skip header
  string dummy;
  input >> dummy;
 
  // determine number of lines of file
  int    n;
  long   lines = 0;
  number x;
  while (!input.eof())
    {
      for(n=0;n<cols;n++)
	input >> x;

      lines+=(!input.eof());
    }
  input.close();

  // set up result matrix
  int newCols = columns+1;
  int newRows = (rows>lines ?rows : lines);
  matrix result(newCols,newRows);
  result.zero();

  int m;
  if (rows&&columns)
    for(n=0;n<rows;n++)
      for(m=0;m<columns;m++)
	result.load(m,n,get(m,n));

  // now add new column from file
  number v=0.;
  ifstream input2;
  input2.open(filename);
  if (input2.fail())
    {
      cerr << "matrix: could not reopen file "
	   << filename
	   << endl;
      
      exit(1);
    }

  // skip header, again...
  input2 >> dummy;

  lines=0;
  while (!input2.eof())
    {
      for(n=1;n<=cols;n++)
	{
	  input2 >> x;
	  if (n==col)
	    v = x;
	}

      if (!input2.eof()) 
	result.load(newCols-1,lines,v);
      lines++;
    }
  input2.close();

  *this=result;
}
/*
number matrix::lndet()
{
  matrix U,W;
  eigenvalvect(*this,U,W);
  
  number result=0.;
  for(int n=0;n<rows;n++)
    result+=log(W.get(n));

  return result;
}*/
/*
void  matrix::polyfit(int order,number ZERO)
{
  // number of template functions, Patrick's dumb version
  int size=0;
  for(int n=0;n<=order;n++)
    for(int m=0;m<=n;m++,size++);

  matrix Y(size);
  matrix F(size,size);
  F.zero();
  Y.zero();
  
  // "scalar product i-th template/data"
  for(number n=0,index=0;n<=order;n++)
    for(number m=0;m<=n;m++,index++)
      for(int x=0;x<columns;x++)
			for(int y=0;y<rows;y++)
				// only matrix values larger than zero are used
				if (get(x,y)!=ZERO)
					Y.add(index,exp(log(get(x,y))+m*log(1.+x)+(n-m)*log(1.+y)));

  // "scalar product i-th template/j-th template"
  for(number n1=0,index1=0;n1<=order;n1++)
    for(number m1=0;m1<=n1;m1++,index1++)
      for(number n2=0,index2=0;n2<=order;n2++)
			for(number m2=0;m2<=n2;m2++,index2++)
				for(int x=0;x<columns;x++)
					for(int y=0;y<rows;y++)
						// only matrix values larger than zero are used
						if (get(x,y)!=ZERO)
							F.add(index1,index2,exp((m1+m2)*log(1.+x)+(n1+n2-m1-m2)*log(1.+y)));

  Y = F.inverse3()*Y;

  // fill the missing elements
  zero();
  for(number n=0,index=0;n<=order;n++)
    for(number m=0;m<=n;m++,index++)
      for(int x=0;x<columns;x++)
			for(int y=0;y<rows;y++)
				add(x,y,exp(m*log(1.+x)+(n-m)*log(1.+y))*Y.get(index));
}*/

void matrix::flipx()
{
  for(int y=0;y<=(int)floor(rows/2-1);y++)
    for(int x=0;x<columns;x++)
    {
    	number tmp = get(x,y);
    	load(x,y,get(x,rows-y-1));
    	load(x,rows-y-1,tmp);
    }
}

void matrix::flipy()
{
  for(int x=0;x<=(int)floor(columns/2-1);x++)
    for(int y=0;y<rows;y++)
      {
	number tmp = get(x,y);
	load(x,y,get(columns-x-1,y));
	load(columns-x-1,y,tmp);
      }
}

matrix matrix::subMatrixKeep(vector<int> elements)
{
	if(columns!=rows)
	{
		clog<<"columns and rows are not equal"<<endl;
		exit(1);
	}
	int subMColumns=elements.size();
	//clog<<"columns="<<input.columns<<"  elements.size()="<<elements.size()<<"  subMColumns="<<subMColumns<<endl;
	matrix subM(subMColumns,subMColumns);

	int m=0;
// 	for(int i=0;i<elements.size();i++)
// 		clog<<"elements["<<i<<"]="<<elements[i]<<endl;
	for(int i=0; i<columns;i++)
	{
		int n=0;
		if(Equal(i,elements))
		{
			for(int j=0;j<columns;j++)
			{
				if(Equal(j,elements))
				{
					subM.load(m,n,get(i,j));
// 					clog<<"m="<<m<<"  n="<<n<<" i="<<i<<" j="<<j<<endl;
					n++;
				}
			}
			m++;
		}
	}
	return subM;
}


matrix subMatrix(matrix input,vector<int> elements)
{
	if(input.columns!=input.rows)
	{
		clog<<"columns and rows are not equal"<<endl;
		exit(1);
	}
	int subMColumns=input.columns-elements.size();
	//clog<<"columns="<<input.columns<<"  elements.size()="<<elements.size()<<"  subMColumns="<<subMColumns<<endl;
	matrix subM(subMColumns,subMColumns);

	int m=0;
// 	for(int i=0;i<elements.size();i++)
// 		clog<<"elements["<<i<<"]="<<elements[i]<<endl;
	for(int i=0; i<input.columns;i++)
	{
		int n=0;
		if(notEqual(i,elements))
		{
			for(int j=0;j<input.columns;j++)
			{
				if(notEqual(j,elements))
				{
					subM.load(m,n,input.get(i,j));
// 					clog<<"m="<<m<<"  n="<<n<<" i="<<i<<" j="<<j<<endl;
					n++;
				}
			}
		m++;
		}
	}
	return subM;
}

bool notEqual(int i,vector<int> elements)
{
	bool result=true;
	for(int j=0;j<elements.size();j++)
		if(i==elements[j])
			result=false;
	return result;
}

bool Equal(int i,vector<int> elements)
{
	bool result=false;
	for(int j=0;j<elements.size();j++)
		if(i==elements[j])
			result=true;
	return result;
}

matrix& matrix::readFromASCII(const char* filename, int cols, int rows)
{
	ifstream fhandler(filename);
	
	alarm(fhandler.fail(),"could not open file",filename);
	if(fhandler.fail())
		exit(1);

	// read in; number of rows may be smaller than rows
	clear();
	vector<matrix> store;
	for(int y=0;y<rows&&!fhandler.fail()&&!fhandler.eof();y++)
	{
		matrix tmpmat(cols);
		for(int x=0;x<cols;x++)
		{	  
			string tmp;
			fhandler >> tmp;
			tmpmat.load(x,atof(tmp.c_str()));
		}

      if (!fhandler.fail()&&!fhandler.eof())
			store.push_back(tmpmat);
	}
	fhandler.close();

	// transfer to matrix array
	resize(cols,store.size());
	for(int y=0;y<store.size();y++)
		for(int x=0;x<cols;x++)
			load(x,y,store[y].get(x));

  return *this;
}

matrix& matrix::readFromASCII_marika(const char* filename)
{
  ifstream fhandler(filename);
  
  alarm(fhandler.fail(),"could not open file",filename);
  if(fhandler.fail())
    exit(1);

  // read in; number of rows may be smaller than rows
  clear();

  string line;
  int counter=0;
  int cols_temp=0,cols=0;
  // vector<stringstream> stream_vec;
  vector< vector<number> > temp_vec_vec;
  while(getline(fhandler,line))
  {
    stringstream stream(line);
    // string str=string(stream.peek());
    // if (str=="#")
    // {
    //   //nColumns++;
    //   string line;
    //   getline(fhandler,line);
    //   clog<<line<<endl;
    // }
    // else
    // {
      //stream_vec.push_back(stream);
      cols_temp=0;
      number temp;
      vector<number> temp_vec;
      while(stream>>temp)
      {
        temp_vec.push_back(temp);
        //clog<<temp<<endl;
        cols_temp++;
      }
      if(cols_temp>0)
      {
        cols=cols_temp;
        temp_vec_vec.push_back(temp_vec);
        counter++;
      }
      else
      {
        clog<<"empty line"<<endl;
      }
      //clog<<"columns="<<cols<<endl;
      //clog<<line_vec[counter]<<endl;
       
    // }
  }
  fhandler.close();
  rows=counter;
  //cols=columns;
  clog<<"rows="<<rows<<endl;
  clog<<"columns="<<cols<<endl;

  // transfer to matrix array
  resize(cols,rows);
  for(int y=0;y<rows;y++)
  {
    for(int x=0;x<cols;x++)
    {
      load(x,y,temp_vec_vec[y][x]);
      //clog<<get(x,y)<<"\t";
    }
    //clog<<endl;
  }

  return *this;
}


matrix euler(number a,number b,number c)
 {
  return D1(a)*D2(b)*D1(c);
 }

 matrix D1(number angle)
 {
  matrix a(3,3);

  number d1[] = {_cos(angle),-_sin(angle),_zero,
                 _sin(angle),_cos(angle),_zero,
                 _zero,_zero,_one};
  
  a.load(d1);

  return a;
 }

 matrix D2(number angle)
 {
  matrix a(3,3);

  number d2[] = {_cos(angle),_zero,+_sin(angle),
                 _zero,_one,_zero,
                 -_sin(angle),_zero,_cos(angle)};

  a.load(d2);

  return a;
 }

 matrix D3(number angle)
 {
  matrix a(3,3);

  number d3[] = {_one,_zero,_zero,
                 _zero,_cos(angle),-_sin(angle),
                 _zero,_sin(angle),_cos(angle)};
  
  a.load(d3);

  return a;
 }

 matrix rot3d(const matrix& a)
 {
  number plugIn[] = {0,-a.elements[2],+a.elements[1],
                     +a.elements[2],0,-a.elements[0],
                     -a.elements[1],a.elements[0],0};

  matrix c(3,3);
  c.load(plugIn);
 
  return c;
 }

// gaussj() from numerical recp. CAUTION! HAS BEEN MODIFIED (ARRAYS)!

#define NRANSI
void gaussj(number **a, int n, number **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol=0,irow=0,j,k,l,ll;
	number big,dum,pivinv,tempr;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef NRANSI 
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */


void polarCoords(number x,number y,number z,
                 number& b, number& l,number& r)
{
 number rho = sqrt(x*x+y*y);

 r   = sqrt(x*x+y*y+z*z);

 if ((rho==0.0) && (z>0))    b = 90.0;
  else 
 if ((rho==0.0) && (z==0.0)) b = 0.0;
  else
 if ((rho==0.0) && (z<0))    b = -90.0;
  else                       b = 180.0/pi*_atan(z/rho);

 if ((x==0.0) && (y==0.0))   l = 0.0;
  else
  {
   number phi = 2.0*180.0/pi*_atan(y/(fabs(x)+rho));

   if ((x>=0.0) && (y>=0.0)) l = phi;
    else
   if ((x>=0.0) && (y<0.0))  l = 360.0+phi;
    else                     l = 180.0-phi;
  }
}
/*
matrix sqrt(const matrix& Y0)
{
  const long    MAXSTEPS=500;
  const number  ACCURACY=1E-6;

  /// based on the algorithmn from Denman&Beavers

  // set start parameter
  matrix Yk = Y0;
  matrix Zk = Y0;
  Zk.neutral();

  number deviation;
  long          k=0;
  int long m,n;
  matrix Ykinv, Zkinv;
  do {
    // make inverts
    Ykinv = Yk.inverse();
    Zkinv = Zk.inverse();

    // next iteration
    Yk.add(Zkinv);Yk.div(2.);
    Zk.add(Ykinv);Zk.div(2.);

    // calculate deviation from desired result
    Ykinv = Yk*Zk;
    deviation=0.;
    for(n=0;n<Y0.rows;n++)
      for(m=0;m<Y0.rows;m++)
	if (m!=n)
	  deviation+=
	    fabs(Ykinv.get(m,n)*Ykinv.get(m,n));
	else
	  deviation+=
	    fabs((Ykinv.get(m,n)-1.)*(Ykinv.get(m,n)-1.));
  }
  while (k++<MAXSTEPS && sqrt(deviation)>ACCURACY);
  
  if (deviation>ACCURACY)
    clog << "sqrt: was not able to achieve desired accuracy "
	 << "[" 
	 << deviation
	 << "]" <<endl;

  return Yk;
}*/
/*
matrix sqrtsym(const matrix & M)
{
  matrix result = M;
  matrix transf = result.diagonal();

  number x;
  int    n;

  for(n=0;n<result.rows;n++)
    if ((x=result.get(n,n))>=0.)
      result.load(n,n,sqrt(x));
    else
      {
	
       //clog << "sqrtsym: WARNING! Negative eigenvalues, setting to zero" << endl;
       //matrix t=M;
       //t.printOut(clog,2);
	

       result.load(n,n,0.);
      }

  return transf*result*transf.t();
}*/
  
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
/*
void jacobi(number **a, int n, number d[], number **v, int *nrot)
{
	int j,iq,ip,i;
	number tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=Vector(1,n);
	z=Vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (number)(fabs(d[ip])+g) == (number)fabs(d[ip])
					&& (number)(fabs(d[iq])+g) == (number)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((number)(fabs(h)+g) == (number)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE

void eigsrt(number d[], number **v, int n)
{
	int k,j,i;
	number p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}


void eigenvalvect(matrix& M, matrix& eigenvectors,matrix& eigenvalues)
{
  if (M.rows!=M.columns)
    nrerror("eigenvalvect: matrix is not square matrix!");

  int N = M.columns;

  eigenvectors.resize(N,N);
  eigenvalues.resize(1,N);

  number** a = Matrix(1,N,1,N);
  number*  d = Vector(1,N);
  number** v = Matrix(1,N,1,N);
  int nrot,m,n;
  for(n=1;n<=N;n++)
    for(m=1;m<=N;m++)
      a[m][n] = M.get(m-1,n-1);

  ///only for symmetric matrices
  jacobi(a,N,d,v,&nrot);
  //eigsrt(d,v,N);

  for(n=1;n<=N;n++)
    {
      for(m=1;m<=N;m++)
	eigenvectors.load(m-1,n-1,v[n][m]);

      eigenvalues.load(n-1,d[n]);
    }

  free_vector(d,1,N);
  free_Matrix(a,1,N,1,N);
  free_Matrix(v,1,N,1,N);
}*/
/*
bool makesvd(matrix& A,matrix& U,matrix& W,matrix& V)
{
  if (!A.size())
    return false;
  
  int N = A.columns;
  int M = A.rows;

  U.resize(N,M);
  W.resize(1,N);
  V.resize(N,N);
  
  number **a = Matrix(1,M,1,N);
  number  *w = Vector(1,N);
  number **v = Matrix(1,N,1,N);
  for(int n=1;n<=N;n++)
    for(int m=1;m<=M;m++)
      a[n][m] = A.get(m-1,n-1);

  // perform the singular value decomposition
  svdcmp(a,M,N,w,v);

  bool zero=false;
  for(int n=1;n<=N;n++)
    {
      for(int m=1;m<=M;m++)
	U.load(m-1,n-1,a[n][m]);

      for(int m=1;m<=N;m++)
	V.load(m-1,n-1,v[n][m]);

      W.load(n-1,w[n]);
      if (w[n]==0.)
	zero=true;
    }

  free_vector(w,1,N);
  free_Matrix(a,1,M,1,N);
  free_Matrix(v,1,N,1,N);

  return zero;
}

/*
#define NRANSI
#define TINY 1.0e-20;

void ludcmp(number **a, int n, int *indx, number *d)
{
	int i,imax=0,j,k;
	number big,dum,sum,temp;
	number *vv;

	vv=Vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}
#undef TINY
#undef NRANSI*/
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */
/*
void lubksb(number **a, int n, int *indx, number b[])
{
	int i,ii=0,ip,j;
	number sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */

/*
#define NRANSI

number pythag(number a, number b)
{
	number absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */
/*
void svdcmp(number **a, int m, int n, number w[], number **v)
{
	int flag,i,its,j,jj,k,l,nm;
	number anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=Vector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((number)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((number)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((number)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}
#undef NRANSI*/
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */
/*
void choldc(number **a, int n, number p[])
{
	void nrerror(char error_text[]);
	int i,j,k;
	number sum;

	for (i=1;i<=n;i++) {
		for (j=i;j<=n;j++) {
			for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					nrerror("choldc failed");
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software ,5+?`[K. */

/*
void choldc_complex(cnumber *a, int n, cnumber *p)
{
	int i,j,k;
	cnumber sum = 0;

	for (i=0;i<n;i++) 
	{
		for (j=i;j<n;j++) 
		{
			for (sum=a[i*n+j],k=i-1;k>=0;k--) sum -= a[i*n+k]*conj(a[j*n+k]);
				if (i == j) 
				{
					alarm(abs(sum)==0,"complex Cholesky failed");
					p[i]=sqrt(sum);
				} 
				else a[j*n+i]=sum/p[i];
		}
	}
}*/
/*
void complexCholesky(matrix& mat_real,matrix& mat_imag)
{
  alarm(mat_real.columns!=mat_imag.columns||mat_real.rows!=mat_imag.rows,
	"real and imaginary part matrices are not of same dimensions in complexCholesky");

  int N = mat_real.columns;
  cnumber* a = new cnumber[N*N];
  cnumber* p = new cnumber[N];

  for(int n=0;n<N;p[n]=0.,n++)
    for(int m=0;m<N;m++)
	a[n+m*N] = cnumber(mat_real.get(n,m),mat_imag.get(n,m));

  choldc_complex(a,N,p);

  mat_real.zero();
  mat_imag.zero();
  for(int n=0;n<N;n++)
    for(int m=n;m<N;m++)
    {
    	mat_real.load(n,m,n==m?real(p[n]):real(a[n+m*N]));
    	mat_imag.load(n,m,n==m?imag(p[n]):imag(a[n+m*N]));
    }

  delete[] a;
  delete[] p;
}

void complexInverse(matrix& mat_real, matrix& mat_imag)
{
  alarm(mat_real.columns!=mat_imag.columns||mat_real.rows!=mat_imag.rows,
	"real and imaginary part matrices are not of same dimensions in complexInverse");

  int N = mat_real.columns;
  
  // load real and imaginary part into one matrix that is inverted
  matrix mat(2*N,2*N);
  for(int n=0;n<N;n++)
    for(int m=0;m<N;m++)
    {
		mat.load(n,m,mat_real.get(n,m));
		mat.load(n+N,m+N,mat_real.get(n,m));
		mat.load(n+N,m,+mat_imag.get(n,m));
		mat.load(n,m+N,-mat_imag.get(n,m));
    }

  mat.invert3();

  // get imaginary part and real part from that matrix
  for(int n=0;n<N;n++)
    for(int m=0;m<N;m++)
    {
    	mat_real.load(n,m,mat.get(n,m));
    	mat_imag.load(n,m,mat.get(n+N,m));
    }
  
  //DONE!
}*/

matrix matrix::ScalarProduct(const matrix& mat2)
{
	if(rows==mat2.rows && columns==mat2.columns)
	{
		matrix product_mat(columns,rows);
		//clog<<"rows="<<rows<<"  columns="<<columns<<endl;
		for(int r=0; r<rows; r++)
			for(int c=0; c<columns; c++)
				product_mat.load(c,r,get(c,r)*mat2.get(c,r));
		return product_mat;
	}
	else
	{
		matrix_error("dimentions don't match for scalar product");
	}
}

matrix matrix::diag()
{
	if(rows!=columns)
		matrix_error("not a square matrix, can't take the diagonal elements");
	matrix vector_mat(rows);
	for(int i=0; i<rows; i++)
		vector_mat.load(i,get(i,i));
	return vector_mat;
}

matrix matrix::power(number power)
{
	matrix power_mat(columns,rows);
	for(int i=0; i<columns; i++)
		for(int j=0; j<rows; j++)
			power_mat.load(i,j,pow(get(i,j),power));
	return power_mat;
}
