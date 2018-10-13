#ifndef FUNCTION_H
#define FUNCTION_H

using namespace std;

#include "globaldef.h"
//#include "nrutil.h"
#include "ticker.h"
#include "matrix.h"
#include "tostring.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


typedef number (WIREFUNC)(number);

const string extension=".table";

/** General handler for functions, implements also data tables (can be stored on disk) and extra/inter-polation

 Allows an easy, relatively easy, way to implement functions. This class includes methods to
 calculate data tables from an arbitrary function that can be stored on disk or restore from
 disk, provided there is a defined interface for a single-argument function.
 The function can then be told to extra/-interpolate (or just one of these) from this data table,
 saving in some situations computation time.

 Moreover, integrals of this function can be calculated, if it, again, can be boiled down to
 a single-argument function.

 Implementations of concrete functions have to be derived from this parent class.

 CAUTION: since 14-01-2005 logarithmic tables are just half-log (only x-axis)
*/

class function
{
 public:

  static const int NONAMECOUNTER   = 1;
  static const int WITHNAMECOUNTER = 2;

    function() 
    {
    	clog<<"in function()"<<endl;
		init();
		clog<<"*******************init"<<endl;
		setName("abstract",WITHNAMECOUNTER);
    }
    function(const char*);
    function(const function&);
    void operator=(const function&);
    virtual ~function();

    void clear();

    string tablename();
	
    void readfromascii(const char*,int,int,int);
    void makeTable(number,number,int,bool=false);
	//makes tables of variable data points for different ranges
    void makeTable_parts(vector<number>,vector<number>,vector<int>,int numberOfParts,bool=false);
    void useTable();
    void noTable();
    bool saveTable();
    bool saveTable_parts();
    bool loadTable(number,number,int,bool=false);
    bool loadTable_noheader(number from, number to, int N,bool logarithmic);
	//loads tables of variable data points
    bool loadTable_parts(vector<number> from, vector<number> to, vector<int> N,int numberOfParts,bool logarithmic=false);
    
    void loadWithValues(matrix&, matrix&,bool=false);
    void loadWithValues(vector<number>&,vector<number>&,bool=false);
    void loadWithValues(matrix&, number, number, bool=false);

    /// use these methods to get values of function
    number value(number);
	//for gsl this is the format
	virtual number integrant(number x, void* params){return value(x);}
	//this was for NR
    virtual number integrant(number x) {return value(x);}
    
    /// alternatives if function depends on more than one variable
    virtual number value(number,number){return 0.0;}
    virtual number value(number,number,number) {return 0.0;}
    virtual number value(number,number,number,number) {return 0.0;}

    /// finds a root of function (one variable) via Newton's method
    number findroot(number x, number dx, number alpha=0.1);
  
    void integrantValues(number,number,int);

    void printOut(number,number,int);
    void printOut(const char*,number,number,int);

    void   setMinMaxReference(number,number);
    number maxReference();
    number minReference();
	 int ShowTableSize();

    void   setName(const char*, int=WITHNAMECOUNTER);
    string getName() {return myName;}

    /// extrapolate over given range for tabelised functions?
    void extrapolationOn();

    /// extrapolate over given range for tabelised functions?
    void extrapolationOff();

    /// take closest value in table as result, no interpolation/extrapolation!
    void interpolationOff();

    /// interpolate (spline)
    void interpolationOn();
    
    /// multiply result with some constant value to rescale
    void multiplyWith(number x) {mult=x;}

    void init();

 protected:

    void copy(const function&);

    int count;

    string   myName;
    bool     noCount;
    bool     takeTable;
    bool     logTable;
    bool     extrapolation;
    bool     interpolation;
    int      tableSize;
    int      tableParts;
    number   tableStart;
    number   tableEnd;
    number   mult;
    number*  referenceY;
    number*  referenceX;
    number*  referenceY2;
	//gsl accelarator object, an itarator for interpolation lokkups
	gsl_interp_accel* acc;
	//gsl spline object
	gsl_spline* spline;

    number interpolate(number);
    number extrapolate(number);
    number find_closest_value(number);

    void   prepareLookupTable();

    /// use this method to define function (overwrite)
    virtual number get(number);
};

/**
   Wires an object of class function to an external function of type 
   (number)(number x) which, for example, could then easily be tabled 
   which then could be used for interpolation later on.
*/
class wiredfunction : public function
{
 public:

  wiredfunction();

  /// wires the object with an external function of type (number)(number)
  void   wirewith(WIREFUNC&);
  number get(number);

 private:
  
  FUNC fct;
};

#endif
