//Global.h
#ifndef __Global_H__
#define __Global_H__
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "RAND\newran.h"
#include <assert.h>
#define PI 3.1415926535898
#define E 2.71828182845904523536
#define Num_Change_Type 6
#define Num_Run 2
#define SWARM_NUM 5
//***********************Enviroment change types in GDBG system ******************************//
enum Change_type{small_step=0, large_step,u_random,recurrent,chaotic,recurrent_noisy};
//********************* Basic component funtion used in composition DBG **********************//
enum Fun_Name{Sphere=0,Rastrigin,Weierstrass,Griewank,Ackley};
//************** flag of minimization or maximization problem used in real space *************//
enum Compare{MIN=0,MAX};
//************** different kind of PSO algorithms *************//
enum swarm_type{GPSO=0,LPSO,APSO,BBPSO,CLPSO};

class General_DBG;

struct Change_Type{
	Change_type type;
	int counter;
};
struct Boundary{						//***************************************************//
	double upper;							//Dimension boudary in the landscape
	double lower;							//**************************************************//
public:
	void Set_Boundary(const double l, const double u){
		lower=l;
		upper=u;
	};
	Boundary & operator=(const Boundary &b){
		if(this==&b) return *this;
		upper=b.upper;
		lower=b.lower;
		return *this;
	};
};
struct Point{
public:
	double x,y;
	double Distance(const Point &p)const{
		return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y));
	}
	Point &operator=(const Point &p){
		if(this==&p) return *this;
		x=p.x;
		y=p.y;
		return *this;
	}
};
class Global{							
public:													//global class : save global variables
	static int num_dim;									// number  of dimensions
	static Cauchy cauchy;								// cauchy random number
	static Normal normal;								// gaussian random number
	static Uniform uniform;								// random number of uniform distribution
	static Change_Type change_type[Num_Change_Type];	// counter for each change style
	static int Swarm_type[SWARM_NUM];
	static int num_peakorfun ;							// number of peaks in Rotation_DBG , number of function in Composition_DBG
	static double hn_s;									// constant number for noralizing all basic function with similar height
	static Boundary boundary;							// search range 
	static double min_height, max_height;				// peak height 
	static double chaotic_constant;						// parameter for chaotic system, between(1,4)
	static double min_width, max_width;					// peak width
	static Change_type change;							// change type
	static int periodicity;								// definite period for values repeating
	static int min_dimension, max_dimension;			// max and min dimension 
	static int change_frequency;						// number of evaluations between two successive changes
	static General_DBG *g_dbg;
	static bool flag_dimension_change;					// flag of dimension change.
	static Compare optimization_type;					// minimization or maximization optimization problem
	static int num_change;								// the number of changes
	static double alpha;								// to control step severity
	static double max_alpha;
	static const int sample_frequency=100;				// frequency of sampling test points during one change
	static int Max_gen;
	static int generation;
	//char resultdir[50];
};

void Initialize_RandomArray(int * a,const int &dim);	// generate a set of radom numbers from 0-|a| without repeat
double Standard_Change(const Change_type,const double min, const double max);
inline double Chaotic_Value(const double x, const double min, const double max){
	// return a value calculated by logistics function 
	if(min>max) return -1;
	double chaotic_value;
	chaotic_value=(x-min)/(max-min);
	chaotic_value=Global::chaotic_constant*chaotic_value*(1-chaotic_value);
	return min+chaotic_value*(max-min);
}
inline double Sin_Value_Noisy(const int x,const double min, const double max, const double amplitude, const double angle,const double noisy_severity=1.){
	// return a value in recurrent with noisy dynamism environment 
	double y;	
	double noisy,t;
	y=min+amplitude*(sin(2*PI*(x+angle)/Global::periodicity)+1)/2.;
	noisy=noisy_severity*Global::normal.Next();
	t=y+noisy;
	if(t>min&&t<max) y=t;
	else y= t-noisy;
	return y;
}

inline double rnorm2(double U,double del)
{
	double u1,u2,x;
	do{u1=Global::uniform.Next();}while(u1==0||u1==1);
	//do{u1=random(0,1);}while(u1==0||u1==1);
	u2=Global::uniform.Next();
	x=sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
	return del*x+U;
}

template <class T>
void Copy(T* destination, const T* source,const int & dim){
	//copy function 
	if(destination==source) return;
	for(int i=0;i<dim;i++)
		destination[i]=source[i];
}
template <class T>							
T Extremum(T * v,const int & size,const Compare & type){
	// return min or max value of set V

	T extreme;
	extreme=v[0];
	//int index=0;
	if(type==MAX){ 
		for(int i=1;i<size;i++){
			if(v[i]>extreme) {
				extreme=v[i];
				//index=i;
			}
		}
	}
	else if(type==MIN){
		for(int i=1;i<size;i++){
			if(v[i]<extreme) {
				extreme=v[i];
				//index=i;
			}
		}
	}
	return extreme;
}

#endif
