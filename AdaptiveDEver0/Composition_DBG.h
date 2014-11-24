//Composition_DBG.h
//#pragma once
#ifndef __Composition_DBG_H__
#define __Composition_DBG_H__
#include "Real_DBG.h"

class Composition_DBG:public Real_DBG{
private:												// *******************class for composition dynamic benchmark generator****************//
	static Composition_DBG* cr_dbg;						// pointer of Composition_DBG class, single instance class
	Boundary * *com_boundary;							// boundary of component functions
	double * converge_severity ;						// severity of converge range for each function
	double * stretch_severity;							// severity of stretching original function, greater than 1 for stretch
														// less than 1 for compress the original function
	double height_normalize_severity;					// constant number for noralizing all basic function with similar height
	Fun_Name *component_function;						// which basic function used to compose the composition function
	
	static const int num_basic_fun=5;					// number of basic functions							
private:
	Composition_DBG(const int num,const int dim);
	void Set_ComBoundary();
public:
	~Composition_DBG();
	static Composition_DBG * Get_Composition_DBG(){
		if(!cr_dbg) cr_dbg= new Composition_DBG(Global::num_peakorfun,Global::num_dim);
		return cr_dbg;
	}
	static void Delete_Composition_DBG(){
		if(cr_dbg){
			delete cr_dbg;
			cr_dbg=0;
		}
	};
	Composition_DBG &operator =(const Composition_DBG &);
	virtual void Set_Rotation_Matrix();					//randomly generate rotation matrx for each basic function
	void Set_Coverge_Sevrity(const double* cs);
	void Set_Stretch_Severity();
	void Set_Basic_Function(const Fun_Name *bf);		//component functions to compose the search space 		
	virtual double Evaluation(const double *const x);
	void Correction(const Fun_Name &f);					// make genes within search range after rotation
	virtual void Calculate_Global_Optima();
	virtual void Parameter_Setting(General_DBG &g_dbg);
public:													// basic five functions
	double F_Sphere();
	double F_Rastrigin();
	double F_Weierstrass();
	double F_Griewank();
	double F_Ackley();
	double Select_Fun(const Fun_Name &f);
public:													// six dynamic changes types
	virtual void Random_Change();
	virtual void Small_Step_Change();
	virtual void Large_Step_Change();
	virtual void Recurrent_Change();
	virtual void Chaotic_Change();
	virtual void Recurrent_Noisy_Change();	
														// dimension changes(linear inrease or decrease)
	virtual void Dimension_Increase();			
	virtual void Dimension_Decrease();
};

#endif

