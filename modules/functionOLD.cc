#include "function.h"

function::function(const function& me)
{
  // this copy constructor DOES NOT call the standard constructor!
  clog<<"In This other constructor"<<endl;
  init();
  copy(me);
}

void function::operator=(const function& me)
{
  // works on an already existing object
  copy(me);
}


void function::init()
{
  //static int counter = 0;
	clog<<"in init"<<endl;
  //clog<<"couner="<<counter<<endl;
  //count = counter++;

  referenceX     = NULL;
  referenceY     = NULL;
  referenceY2    = NULL;
  takeTable      = false;
  tableStart     = 0.0;
  tableEnd       = 0.0;
  tableSize      = 0;
  logTable       = false;
  extrapolation  = true;
  interpolation  = true;
  mult           = 1.;
  clog<<"initialized"<<endl;
}

void function::setName(const char* name,int param)
{
  myName = name;
  clog<<"in setname, my name is: "<<myName<<endl;
  switch (param) 
  {
    case WITHNAMECOUNTER:
      noCount = false;
      break;
    case NONAMECOUNTER:
      noCount = true;
      break;
    default:
      noCount=false;
      break;
  }
}
       
function::function(const char* name)
{
  init();  
  myName  = name;
  noCount = true;
}

function::~function()
{
  clear();
}

void function::clear()
{
  if (referenceX!=NULL)
    delete[] referenceX;
  if (referenceY!=NULL)
    delete[] referenceY;
  if (referenceY2!=NULL)
    delete[] referenceY2;

  init();
}

void function::useTable()
{
    if (referenceX!=0 && referenceY!=0)
	takeTable = true;
}

void function::noTable()
{
    takeTable = false;
}

string function::tablename()
{
  return noCount ?
    string(tablePath)+myName+extension :
    string(tablePath)+myName+string("_")+toString(count)+extension;
}

number function::maxReference()
{
  if (!takeTable)
    return 0.;

  // return maximum x-value of reference table
  return logTable ? exp(tableEnd) : tableEnd;
}

int function::ShowTableSize()
{
  if (!takeTable)
    return 0.;

  // return maximum x-value of reference table
  return tableSize;
}


number function::minReference()
{
  if (!takeTable)
    return 0.;

  // return minimum x-value of reference table
  return logTable ? exp(tableStart) : tableStart;
}

void function::setMinMaxReference(number min,number max)
{
  if (!takeTable)
    return;

  tableStart = (logTable ? log(min) : min);
  tableEnd   = (logTable ? log(max) : max);
}

/*
   -----------------------------------------------
   ---- MAKE DATA TABLE --------------------------
   -----------------------------------------------
*/
void function::makeTable(number from, number to, int N, bool logarithmic)
{
	delete[] referenceX;
	delete[] referenceY;

	referenceX	= new number[N];
	referenceY	= new number[N];
	tableStart	= from;
	tableEnd		= to;
	tableSize		= N;
	takeTable		= true;

	if (logarithmic&&(from>0))
	{
		logTable   = true;
		tableStart = log(tableStart);
		tableEnd   = log(tableEnd);
	}
	else 
	logTable = false;

	int n;
	number   x;

	ticker clock;
	for(n=0;n<N;n++)
	{
		if (n>0 && !(n%((int) (N/20))))
			clock.tick(1.*n/N);

		x = tableStart+(tableEnd-tableStart)/(N-1.0)*n;
		if (!logTable)
		{
			referenceY[n] = get(x);
			referenceX[n] = x;
			
		}
		else
		{
			referenceY[n] = get(exp(x));
			referenceX[n] = x;
		}
	}

	prepareLookupTable();

	clog << "\r";
	for(int i=0;i<50;i++)
	clog << " ";
	clog << endl;
}

void function::makeTable_parts(vector<number> from,vector<number> to, vector<int> N,int numberOfParts, bool logarithmic)
{
    delete[] referenceX;
    delete[] referenceY;

    tableSize=0;
    for(int i=0;i<numberOfParts;i++)
	   tableSize+= N[i];

    referenceX  = new number[tableSize];
    referenceY  =  new number[tableSize];
    tableStart    = from[0];
    tableEnd      = to[numberOfParts-1];
//     tableSize      = N[i];
    takeTable      = true;

    if (logarithmic&&(tableStart>0))
    {
    	logTable   = true;
    	tableStart = log(tableStart);
    	tableEnd   = log(tableEnd);
    	for(int i=0; i<numberOfParts; i++)
    	{
    		from[i]=log(from[i]);
    		to[i]=log(to[i]);
    	}
    }
    else 
      logTable = false;

    number   x;
    int Nsum=0;
    ticker clock;
    for(int i=0; i<numberOfParts;i++)
    {
    	//cout<<"#i="<<i<<endl;
    	for(int n=0;n<N[i];n++)
    	{
    		int point=n+Nsum;
    		if (n>0 && !(n%((int) (tableSize/20))))
    		  clock.tick(1.*n/tableSize);

    		x = from[i]+(to[i]-from[i])/(N[i]-1.0)*n;

    		if (!logTable)
    		{
    		    referenceY[point] = get(x);
    		    referenceX[point] = x;
    		}
    		else
    		{
    			//referenceY[n+1] = log(get(exp(x)));
    			referenceY[point] = get(exp(x));
    			referenceX[point] = x;
    			//cout<<exp(x)<<endl;
    		}
	   }
	   Nsum+=N[i];
  }

  prepareLookupTable();

//     clog << "\r";
//     for(int i=0;i<50;i++)
//       clog << " ";
//     clog << endl;
}

/*
   -----------------------------------------------
   ---- TAKE FUNCTION VALUE ----------------------
   -----------------------------------------------
*/
number function::value(number x)
{
  static number y;
  //--- decide whether to calculate value or to interpolate from reference table

  // no table -> evaluate value according to function's rule
  if (!takeTable)
    y= mult*get(x);//mult=1;
  // reference table, but x outside interpolation range ?
//   else if ((!logTable && ((x<tableStart) || (x>tableEnd))) ||
// 	   (logTable && ((log(x)<tableStart) || (log(x)>tableEnd))))
	else if ((!logTable && ((x<tableStart) || (x>tableEnd))) ||
	   (logTable && ((x<exp(tableStart)) || (x>exp(tableEnd)))))
  {
    y= mult*(logTable ? extrapolate(log(x)) : extrapolate(x));
	 //clog<<"TableStart="<<tableStart<<" TableEnd="<<tableEnd<<endl;
	 //clog<<"in value function, log(x)="<<log(x)<<"  y="<<y<<endl;
	 //clog<<"in value function, log(x)="<<log(x)<<"  y="<<y<<endl;
  }
    //y= mult*(logTable ? exp(extrapolate(log(x))) : extrapolate(x));
  // use reference table via interpolation
  else 
    //y=mult*(logTable ? exp(interpolate(log(x))) : interpolate(x));
    y=mult*(logTable ? interpolate(log(x)) : interpolate(x));
  //clog<<"  y="<<y<<endl;

  return (!isnan(y) ? y : 0.);
}


/*
   -----------------------------------------------
   ---- EXTRAPOLATION ----------------------------
   -----------------------------------------------
*/

void function::extrapolationOn()
{
  extrapolation = true;
}

void function::extrapolationOff()
{
  extrapolation = false;
}

void function::interpolationOff()
{
  interpolation = false;
}

void function::interpolationOn()
{
  interpolation = true;
}


number function::extrapolate(number x)
{
  static number x1,y1;
  static number x2,y2;
  static number x3,y3;

  // handles so far only tables greater than two elements
  
  if (tableSize<3 || !extrapolation)
    return 0.0;

  if (x>referenceX[tableSize])
    {
      x1 = referenceX[tableSize-3];
      y1 = referenceY[tableSize-3];
      x2 = referenceX[tableSize-2];
      y2 = referenceY[tableSize-2];
      x3 = referenceX[tableSize-1];
      y3 = referenceY[tableSize-1];
    }
  else
    {
      x1 = referenceX[0];
      y1 = referenceY[0];
      x2 = referenceX[1];
      y2 = referenceY[1];
      x3 = referenceX[2];
      y3 = referenceY[2];
    }

  return (x-x1)*(x-x2)/(x3-x1)/(x3-x2)*y3 +
         (x-x1)*(x-x3)/(x2-x1)/(x2-x3)*y2 +
         (x-x2)*(x-x3)/(x1-x2)/(x1-x3)*y1;
}


/*
   -----------------------------------------------
   ---- INTERPOLATION ----------------------------
   -----------------------------------------------
*/


number function::interpolate(number x)
{
  static number y;

  // do splint-interpolation (if interpolation is done!)
  if (interpolation)
  {
	  // change this
	   y=gsl_spline_eval(spline,x,acc);
  }
  else
    // just find closest value in table
    y = find_closest_value(x);

  return y;
}


number function::find_closest_value(number x)
{
  number ybest=0.;
  number dbest=0.;
  for(int i=0;i<tableSize;i++)
    if (i==0||fabs(x-referenceX[i])<dbest)
    {
		  dbest = fabs(x-referenceX[i]);
		  ybest = referenceY[i];
    }

  return ybest;
}


number function::get(number x)
{
  // dummy function 

  return 0.0;
}

bool function::saveTable()
{
    if (!tableSize)
	return false;

    /*
    clog << "\r"
	 << tablename()
	 << ": saving reference table to disk " 
	 << endl;
    */

    ofstream file;
    file.open(tablename().c_str());
	 file.precision(20);
    

    file << "# " 
	  << tableStart << " " 
	  << tableEnd   << " " 
	  << tableSize  << " " 
	  << logTable << endl;

    file.width(20);

    int n;
    for(n=0;n<tableSize;n++)
      file << referenceX[n]
	       << "\t"
	       << referenceY[n] 
	       << endl;

    file.close();

    return true;
}

bool function::saveTable_parts()
{
    if (!tableSize)
	return false;

    /*
    clog << "\r"
	 << tablename()
	 << ": saving reference table to disk " 
	 << endl;
    */

    ofstream file;
    file.open(tablename().c_str());

    file << "# " 
      	 << tableStart << " " 
      	 << tableEnd   << " " 
      	 << tableSize  << " " 
      	 << tableParts << " "
      	 << logTable << endl;

    file.precision(20);
    file.width(20);

    int n;
    for(n=0;n<tableSize;n++)
      file << referenceX[n]
      	   << "\t"
      	   << referenceY[n] 
      	   << endl;

    file.close();

    return true;
}

bool function::loadTable(number from, number to, int N,bool logarithmic)
{
	ifstream file;
	file.open(tablename().c_str());
	if (!file.is_open())
	{
		makeTable(from,to,N,logarithmic);
		saveTable();
		return true;
	}

	//clog<<tablename()<< ":restoring reference table from disk ";
	
	delete[] referenceY;
	delete[] referenceX;

	char dummy;

	file.get(dummy);
	file >> tableStart; 
	file >> tableEnd;
	file >> tableSize;
	file >> logTable;
	//clog<<"What's going on????????????????????????"<<endl;

	//clog << " (";

// 	if (logTable)
// 		clog << "log ";
// 	clog<<"[" << tableStart<<","<<tableEnd<<"] ";
// 	clog<< "N=" << tableSize<<")"<<endl;

	referenceY = new number[tableSize];
	referenceX = new number[tableSize];

	int n;
	for(n=0;n<tableSize;n++)
	{
		file >> referenceX[n];
		file >> referenceY[n];
	}

	file.close();
	//clog<<10<<"\t"<<referenceX[10]<<endl;
	//clog<<10<<"\t"<<referenceY[10]<<endl;
	takeTable = true;
	prepareLookupTable();

	return true;
}
///another loadTable which doesn't have a header
bool function::loadTable_noheader(number from, number to, int N,bool logarithmic)
{
    ifstream file;
    file.open(tablename().c_str());
    if (!file.is_open())
      {
	makeTable(from,to,N,logarithmic);
	saveTable();

	return true;
      }

    clog << tablename()
	 << ":restoring reference table from disk ";
	
    delete[] referenceY;
    delete[] referenceX;

//     char dummy;
// 
//     file.get(dummy);
    tableStart=from; 
    tableEnd=to;
    tableSize=N;
    logTable=logarithmic;


    
    if (logTable)
      clog << "log table"<< endl;

    referenceY = new number[tableSize];
    referenceX = new number[tableSize];

    int n;
    for(n=0;n<tableSize;n++)
    {
		  file >> referenceX[n];
		  file >> referenceY[n];
    }

    file.close();

    takeTable = true;

    prepareLookupTable();

    return true;
}

///my loadTable, to make tables with different binning for different arguments
bool function::loadTable_parts(vector<number> from, vector<number> to, vector<int> N,int numberOfParts,bool logarithmic)
{
    ifstream file;
    file.open(tablename().c_str());
    if (!file.is_open())
    {
    	clog<<"writing table:"<<tablename()<<endl;
    	makeTable_parts(from,to,N,numberOfParts,logarithmic);
    	saveTable_parts();
    	return true;
    }

    clog << tablename()<< ":restoring reference table from disk ";
	
    delete[] referenceY;
    delete[] referenceX;

    char dummy;

    file.get(dummy);
    file >> tableStart;
    file >> tableEnd;
    file >> tableSize;
    file >> tableParts;
    file >> logTable;

    clog << " (";

    if (logTable)
      clog << "log ";
    clog << "["  << tableStart 
	       << ","  << tableEnd <<"] ";
    clog << "N=" << tableSize 
	       <<", number of parts="<<tableParts
	       << ")"  << endl;

    referenceY = new number[tableSize];
    referenceX = new number[tableSize];

    int n;
    for(n=0;(n<tableSize)&&(!file.eof());n++)
    {
    	file >> referenceX[n];
    	file >> referenceY[n];
    }

    file.close();

    takeTable = true;

    prepareLookupTable();

    return true;
}

void function::integrantValues(number from,number to,int N)
{
  int n;
  number   x;

  for(n=0;n<N;n++)
    {
      x = from+(to-from)/(N-1.0)*n;
      cout << x << "\t" << integrant(x) << endl;
    }
}

void function::loadWithValues(matrix& x, matrix& y,bool logarithmic)
{
  int size;

  x.columns==1 ? size = x.rows : size = x.columns;
  if (y.columns!=size && y.rows!=size)//??
  {
    cerr << "function: could not load reference table, not equally sized" << endl;
    exit(1);
  }

  delete[] referenceX;
  delete[] referenceY;

  referenceX = new number[size];
  referenceY = new number[size];

  for(int n=0;n<size;n++)
  {
    referenceX[n] = x.elements[n];
    referenceY[n] = y.elements[n];
  }

  takeTable  = true;
  logTable   = logarithmic;
  tableStart = x.min();
  tableEnd   = x.max();
  tableSize  = size;

  prepareLookupTable();
}

void function::loadWithValues(vector<number>& x, vector<number>& y,bool logarithmic)
{
  int size=x.size();
	clog<<"in loadWithValues in function"<<endl;
	clog<<"size is:"<<size<<endl;
  //x.columns==1 ? size = x.rows : size = x.columns;
  if (y.size()!=size)//??
  {
    cerr << "function: could not load reference table, not equally sized" << endl;
    exit(1);
  }

  delete[] referenceX;
  delete[] referenceY;

  referenceX = new number[size];
  referenceY = new number[size];

  for(int n=0; n<size; n++)
  {
    referenceX[n] = x[n];
    referenceY[n] = y[n];
  }
   
  takeTable  = true;
  logTable   = logarithmic;
  tableStart = *min_element(x.begin(), x.end());
  tableEnd   = *max_element(x.begin(), x.end());
  tableSize  = size;

  prepareLookupTable();
}


void function::loadWithValues(matrix& y, number x0, number dx, bool logarithmic)
{
  int size;

  y.columns==1 ? size = y.rows : size = y.columns;

  delete[] referenceX;
  delete[] referenceY;

  referenceX = new number[size];
  referenceY = new number[size];

  number x = x0;
  for(int n=0;n<size;n++)
  {
    referenceX[n] = x;
    referenceY[n] = y.elements[n];
    x+=dx;
  }

  takeTable  = true;
  logTable   = logarithmic;
  tableStart = x0;
  tableEnd   = x-dx;
  tableSize  = size;

  prepareLookupTable();
}

void function::printOut(number x1,number x2,int N)
{
  number dx = (x2-x1)/N;
  number x;

  for(x=x1;x<x2;x+=dx)
    cout << x << "\t" << value(x) << endl;
}

void function::printOut(const char* filename, number x1,number x2,int N)
{
  ofstream filehandler;
  filehandler.open(filename);
  if (filehandler.fail())
  {
    cerr << "function: could not open file " << filename << " -> no data table" << endl;
    return;
  }

  number dx = (x2-x1)/N;
  number x;

  ticker clock;
  int long count = 0;
  for(x=x1;x<x2;x+=dx)
  {
      count+=1;
      if (!fmod(count,floor(N/20.)))
	       clock.tick(1.*count/N);

      filehandler << x << "\t" << value(x) << endl;
  }

  filehandler.close();
}

void function::prepareLookupTable()
{
  // is called once any time a new lookup table is defined
  // makes preparations for Spline algorithmn

  // Spline needs first derivatives at the left and right edge
  // here: approximated by outer 2 points
  clog<<"in prepareLookupTable"<<endl;
  //number dy1 = (referenceY[1]-referenceY[0])/(referenceX[1]-referenceX[0]);
  //clog<<dy1<<endl;
  //number dy2 = (referenceY[tableSize-1]-referenceY[tableSize-2])
//					/(referenceX[tableSize-1]-referenceX[tableSize-2]);
  //clog<<dy2<<endl;
  acc=gsl_interp_accel_alloc();
  //clog<<"accelator allocated"<<endl;
  spline=gsl_spline_alloc(gsl_interp_cspline,tableSize);
  //clog<<"spline allocated"<<endl;
  //clog<<"tableSize="<<tableSize<<"  referenceX[0]="<<referenceX[0]<<"   referenceY[0]="<<referenceY[0]<<endl;
  //number X=referenceX
  gsl_spline_init(spline,referenceX,referenceY,tableSize);
  //clog<<"spline initialized"<<endl;
}

void function::copy(const function& me)
{
  clear();

  takeTable     = me.takeTable;
  tableStart    = me.tableStart;
  tableEnd      = me.tableEnd;
  tableSize     = me.tableSize;
  logTable      = me.logTable;
  noCount       = me.noCount;
  myName        = me.myName;
  extrapolation = me.extrapolation;
  interpolation = me.interpolation;
  mult          = me.mult;

  if (me.referenceX!=NULL)
  {
    referenceX = new number[tableSize];
    referenceY = new number[tableSize];
    referenceY2= new number[tableSize];

    for(int i=0;i<tableSize;i++)
  	{
  	  referenceX[i] = me.referenceX[i];
  	  referenceY[i] = me.referenceY[i];
  	  referenceY2[i]= me.referenceY2[i];
  	}
  }
  else
  {
    referenceX  = NULL;
    referenceY  = NULL;
    referenceY2 = NULL;
  }
}

void function::readfromascii(const char* file,int cols,int colx,int coly)
{
	vector<number> x,y;
	ifstream fhandler;
  
	fhandler.open(file);
	alarm(fhandler.fail()||fhandler.eof(), "function: could not read from ASCII file",file);

	number* tmp = new number[cols];
	while (!fhandler.eof())
	{
		for(int n=0;n<cols;n++)
		{
			string tmp2;
			fhandler >> tmp2;
			tmp[n] = atof(tmp2.c_str());
		}
      
      if (fhandler.eof()||fhandler.fail())
			break;

		x.push_back(tmp[colx]);
		y.push_back(tmp[coly]);
	}
	fhandler.close();
	delete[] tmp;
  
	matrix vecx(x.size());
	matrix vecy(y.size());
	for(int n=0;n<x.size();n++)
  {
		vecx.load(n,x[n]);
		vecy.load(n,y[n]);
  }
	loadWithValues(vecx,vecy);
	extrapolationOff();
}

number function::findroot(number x, number dx, number alpha)
{
  number df;
  int iterations = 0;
  number xold    = x;
  do{
    // derivative
    df = 
      (-value(x+2.*dx)+8.*value(x+dx)-8.*value(x-dx)+value(x-2.*dx))/12./dx;

    xold = x;
    x   -= value(xold)/df*alpha;
  } while (iterations++<1E4&&fabs(x-xold)>1E-5);

  return x;
}

wiredfunction::wiredfunction()
{
  fct = NULL;
}

void wiredfunction::wirewith(WIREFUNC& f)
{
  fct = &f;
}

number wiredfunction::get(number x)
{
  return (*fct)(x);
}
