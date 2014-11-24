//Real_DBG.h
#pragma once
#include "General_DBG.h"
#include "Matrix.h"
#include <fstream>
using namespace std;

class Real_DBG: public General_DBG{
protected:									//*************base class for continuous problems***************//
	double *genes;							// real coded solution
	int num_peakorfun;						// number of peaks in Rotation_DBG , number of function in Composition_DBG 
	double **position;						// positions of local or global optima(local optima in Rotation_DBG,
											// global optima of basic function in Composition_DBG)
	double **initial_position;				// save the initial positions 
	double *height;							// peak height in Rotation_DBG, height of global optima in Composition_DBG
	Boundary *boundary;						// solution space 
	const double min_height;				// minimum height of all peaks(local optima) in Rotation_DBG(Composition_DBG) 
	const double max_height;				// maximum height of all peaks(local optima) in Rotation_DBG(Composition_DBG) 
	double height_severity;
	bool prediction;						// the next change of function can be predictable or not
	double *fit;							// objective value of each basic funciton in Composition_DBG, peak height in Rotation_DBG
	double *weight;							// weight value of each basic function in Composition_DBG,  peak width in Rotation_DBG
	Matrix *rotation_matrix;				// orthogonal rotation matrixes for each function
	double global_optima;					// global optima value
	int ***rotation_planes;					// save the planes rotated during one periodicity
	double *global_optima_position;			// position of global optima
public:
	Real_DBG(const int num_peakorfun,const int dim);
	virtual ~Real_DBG()=0;
	virtual void Position_Change(){};
	virtual void Height_Change(){};
	Real_DBG & operator=(const Real_DBG &r_dbg);
	void Set_Boundary(const Boundary * b);
	void Set_Boundary(const double * b);
	void Set_Boundary(const double &b);
	void Set_Boundary(const Boundary &b);

	double *Get_Genes()const{
		return genes;
	}
	const Boundary * Get_Boundary()const {
		return boundary;
	}
	void Set_Height(const double *h){
		Copy(height,h,num_peakorfun);
	}
	void Set_Position(const double **p){
		for(int i=0;i<num_peakorfun;i++){
			Copy(position[i],p[i],dimension);
			Copy(initial_position[i],p[i],dimension);
		}
	}
	double *Get_Height(){
		return height;
	}
	double ** Get_Position(){
		return position;
	}
	double *Get_Global_Optima_Position(){
		return global_optima_position;
	}
	virtual void Set_Weight(const double w){
		for(int i=0;i<num_peakorfun;i++)
			weight[i]=w;
	}
	void Set_Prediction(const bool p){
		prediction=p;
	}
	virtual void Set_Rotation_Matrix(){};	
	void Set_Height_Severity(const double hs);
	virtual bool Set_Periodicity(const int p);
	virtual double Evaluation(const double *const x){
		return 0.;
	};
	void Correction();
	double Get_Global_Optima();
	virtual void Calculate_Global_Optima(){};
	virtual void Parameter_Setting(General_DBG &g_dbg);
public:								// six dynamic changes types
	void Height_Standard_Change();
	void Position_Standard_Change( double );
	virtual void Random_Change(){};
	virtual void Small_Step_Change(){};
	virtual void Large_Step_Change(){};
	virtual void Recurrent_Change(){};
	virtual void Chaotic_Change(){};
	virtual void Recurrent_Noisy_Change(){};
	virtual void Dimension_Increase(){};
	virtual void Dimension_Decrease(){};
	
public:													// used for debug
	void Print_Fun(ofstream & f){
		for(int i=0;i<1;i++){//num_peakorfun
			f<<"f NO. "<<i<<": "<<height[i]<<" "<<weight[i]<<endl;
			for(int j=0;j<dimension;j++)
				f<<position[i][j]<<" ";
			f<<endl<<endl;
		}
	}
};

