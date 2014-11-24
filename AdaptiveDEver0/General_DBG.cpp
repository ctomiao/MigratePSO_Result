#include "stdafx.h"
#include "General_DBG.h"
#include <iostream>

using namespace std;

//constructor of General_DBG 
General_DBG::General_DBG(const int dim):dimension(dim){
	periodicity=0;				//set default value=0
	changes_counter=0;
	flag_dimension_change=false;
	dir_dimension_change=true;
}
General_DBG::~General_DBG(){

}
General_DBG::General_DBG(const General_DBG  &g){

	dimension=g.dimension;
	changes_counter=g.changes_counter;
	Set_Change_Frequency(g.change_frequency);
	Set_Change_Type(g.change_type);

	Set_Periodicity(g.periodicity);
	
}
 //assignment operator 
General_DBG& General_DBG::operator =(const General_DBG  & g){
	
	if(this==&g) return *this;
	if(dimension!=g.dimension){
		cout<<"dimension should be the same"<<endl;
		exit(0);
	} 
	if(change_type!=g.change_type){
		cout<<"change type should be the same"<<endl;
		exit(0);
	} 
	changes_counter=g.changes_counter;
	Set_Change_Frequency(g.change_frequency);
	Set_Change_Type(g.change_type);
	Set_Periodicity(g.periodicity);
	return *this;
}

const int General_DBG::Get_Dimension()const{
	return dimension;
}
const int General_DBG::Get_Changes_Counter()const{
	return changes_counter;
}

void General_DBG::Set_Change_Frequency(const int f){
	change_frequency=f;
}
const int General_DBG::Get_Change_Frequency() const{
	return change_frequency;
}
void General_DBG::Set_Change_Type(const Change_type t){
	change_type=t;
}
const Change_type General_DBG::Get_Change_Type() const{
	return change_type;
}
bool General_DBG::Set_Periodicity(const int p){
	if(p<1) return false;
	periodicity=p;
	return true;
}
void General_DBG::Set_RR_Severity(const float p){
		recurrent_noisy_servity=p;
}
const int General_DBG::Get_Periodicity()const{
		return periodicity;
}
void General_DBG::Change(){
	changes_counter++;
	switch(Get_Change_Type()){
	case u_random:
		Random_Change();
		break;
	case recurrent:
		Recurrent_Change();
		break;
	case recurrent_noisy:
		Recurrent_Noisy_Change();
		break;
	case small_step:
		Small_Step_Change();
		break;
	case large_step:
		Large_Step_Change();
		break;
	case chaotic:
		Chaotic_Change();
		break;
	default :
		break;
	} 
	
}
void General_DBG::Dimension_Change(){
	if(flag_dimension_change==false) return;

	Global::num_dim=dimension;
	if(Global::num_dim==Global::min_dimension)
		dir_dimension_change=true;
	if(Global::num_dim==Global::max_dimension )
		dir_dimension_change=false;

	if(dir_dimension_change==true) Global::num_dim++;
	else 	Global::num_dim--;
	
	if(Global::num_dim<dimension)	Dimension_Decrease();
	else Dimension_Increase();

}
void General_DBG::Parameter_Setting(General_DBG &g_dbg){
			
	change_frequency=g_dbg.change_frequency;				
	change_type=g_dbg.change_type;			
	periodicity=g_dbg.periodicity;																			
	changes_counter=g_dbg.changes_counter;				
	recurrent_noisy_servity=g_dbg.recurrent_noisy_servity;		
	flag_dimension_change=g_dbg.flag_dimension_change;			
	dir_dimension_change=g_dbg.dir_dimension_change;			
}

