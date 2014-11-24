//Global.cpp
#include "stdafx.h"
#include "Global.h"

using namespace std;

int Global::num_dim=10;
Cauchy Global::cauchy;
Normal Global::normal;
Uniform Global::uniform;
Change_Type Global::change_type[Num_Change_Type];
int Global::num_peakorfun;
double Global::hn_s;
Boundary Global::boundary;
double Global::min_height;
double Global::max_height;
double Global::chaotic_constant;
double Global::min_width;
double Global::max_width;
Change_type Global::change;
int Global::periodicity;
int Global::min_dimension;
int Global::max_dimension;
General_DBG *Global::g_dbg;
bool Global::flag_dimension_change=false;
Compare Global::optimization_type;
double Global::alpha=0.04;
double Global::max_alpha=0.1;
int Global::change_frequency=10000;
//int Global::num_change=60;// the original
int Global::num_change=20;
int Global::Max_gen;
int Global::generation;
void Initialize_RandomArray(int * a,const int &dim){
	int * temp=new int[dim];
	for(int i=0;i<dim;i++)	temp[i]=i;
	int d=dim;
	for(int i=0;i<dim;i++){
		int t= (int)(d*Global::uniform.Next());
		a[i]=temp[t];
		for(int k=t;k<d-1;k++)
			temp[k]=temp[k+1];
		d--;
	}
	delete []temp;
}
double Standard_Change(const Change_type T, const double min, const double max){
	double step,sign;
	switch(T){
	case small_step:
		step=-1+2*Global::uniform.Next();
		step=Global::alpha*step*(max-min);
		break;
	case u_random:
		step=Global::normal.Next();
	case large_step:
		step=-1+2*Global::uniform.Next();
		if(step>0)sign=1;
		else if(step<0) sign=-1;
		else sign=0; 
		step=(Global::alpha*sign+(Global::max_alpha-Global::alpha)*step)*(max-min);
		break;
	case recurrent:
	case chaotic:
	case recurrent_noisy:
		break;
	}
	return step;
}