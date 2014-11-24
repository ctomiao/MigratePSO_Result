// General_DBG.h
#pragma once

#include "Global.h"

class General_DBG{						//****************************************************//
protected:								// basic class for general dynamic benchmark gennerator		
	//double pheno_severity;				// change probability of dimension pheno_severity is (0,1)
	int change_frequency;				// number of evaluations between two successive changes
	Change_type change_type;			// type of the current change
	int periodicity;					// definite period for values repeating									
	int dimension;						// problem dimension
	int changes_counter;				// counter of number of changes
	float recurrent_noisy_servity;		// deviation servity from the trajactory of recurrent change
	bool flag_dimension_change;			// flag=true, the number of dimensions change, otherwise no change, 
										// default value is false
	bool dir_dimension_change;			// direction of change, dir=true means increasing the dimension, otherwise decrease it
public:
	General_DBG(const int dim);
	virtual ~General_DBG()=0;
	General_DBG(const General_DBG  &);
	General_DBG & operator=(const General_DBG  &);

	const int Get_Dimension()const;
	const int Get_Changes_Counter()const;
	void Set_Change_Frequency(const int f);
	const int Get_Change_Frequency() const;
	void Set_Change_Type(const Change_type t);
	const Change_type Get_Change_Type() const;
	virtual bool Set_Periodicity(const int p);
	void Set_RR_Severity(const float p);
	const int Get_Periodicity()const;
	virtual void Parameter_Setting(General_DBG &g_dbg);
	void Set_Dimension_Change(const bool f){
		flag_dimension_change=f;
	}
	const bool Get_Dimension_Change_Flag(){
		return flag_dimension_change;
	}
	const bool Get_Dimension_Change_Dir(){
		return dir_dimension_change;
	}
public:											// six dynamic changes types
	void Change();
	virtual void Random_Change(){};
	virtual void Small_Step_Change(){};
	virtual void Large_Step_Change(){};
	virtual void Recurrent_Change(){};
	virtual void Chaotic_Change(){};
	virtual void Recurrent_Noisy_Change(){};
												// dimension changes(linear inrease or decrease)
	void Dimension_Change();
	virtual void Dimension_Increase(){};
	virtual void Dimension_Decrease(){};
};

