//Rotation_DBG.h
// Programmed 
#pragma once
# include "Real_DBG.h"
class Rotation_DBG:public Real_DBG{						//***************Rotation Dynamic Benchmark ******************//
private:												//
	static Rotation_DBG * rot_dbg;						// pointer of Rotation_DBG class, single instance class
	double width_severity;								// width severity of each peak
	const double min_width, max_width;					// maximun and minimum width for all peaks 
private:
	Rotation_DBG(const int num,const int dim);
public:
	~Rotation_DBG();
	static Rotation_DBG* Get_Rotation_DBG(){
		if(!rot_dbg) rot_dbg=new Rotation_DBG(Global::num_peakorfun,Global::num_dim);
		return rot_dbg;
	}
	static void Delete_Rotation_DBG(){
		if(rot_dbg){
			delete rot_dbg;
			rot_dbg=0;
		}
	}
	Rotation_DBG& operator=(const Rotation_DBG &);
	virtual void  Set_Weight(const double w){
		for(int i=0;i<num_peakorfun;i++)
			if(Get_Change_Type()==chaotic)
				weight[i]=min_width+(max_width-min_width)*Global::uniform.Next();
			else
				weight[i]=w;
	};
	void Set_Width_Severity(const double sw);
	virtual bool Set_Periodicity(const int p);
	virtual double Evaluation(const double *const x);
	virtual void Calculate_Global_Optima();
	virtual void Parameter_Setting(General_DBG &g_dbg);
public:													// six dynamic changes types
	void Width_Standard_Change();
	virtual void Random_Change();
	virtual void Small_Step_Change();
	virtual void Large_Step_Change();
	virtual void Recurrent_Change();
	virtual void Chaotic_Change();
	virtual void Recurrent_Noisy_Change();
	virtual void Dimension_Increase();
	virtual void Dimension_Decrease();
};

