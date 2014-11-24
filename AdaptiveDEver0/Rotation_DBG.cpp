//Rotation_DBG.cpp
#include "stdafx.h"
#include "Rotation_DBG.h"
#include "Global.h"
Rotation_DBG* Rotation_DBG::rot_dbg=0;

Rotation_DBG::Rotation_DBG(const int num, const int dim):Real_DBG(num,dim),min_width(Global::min_width),max_width(Global::max_width){
	width_severity=0;
}
Rotation_DBG::~Rotation_DBG(){
	
}
Rotation_DBG& Rotation_DBG::operator =(const Rotation_DBG & r_dbg){
	if(this==& r_dbg) return *this;
	Real_DBG::operator =(r_dbg);
	width_severity=r_dbg.width_severity;
	return *this;
}
void Rotation_DBG::Random_Change(){
	Height_Standard_Change();
	Width_Standard_Change();
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[u_random].counter++;
}
void Rotation_DBG::Recurrent_Change(){
	
	double initial_angle;
	double height_range=max_height-min_height;
	double width_range=max_width-min_width;
	for(int i=0;i<num_peakorfun;i++){
		initial_angle=(double)Get_Periodicity()*i/num_peakorfun;
		height[i]=min_height+height_range*(sin(2*PI*(Global::change_type[recurrent].counter+initial_angle)/Get_Periodicity())+1)/2.;
		weight[i]=min_width+width_range*(sin(2*PI*(Global::change_type[recurrent].counter+initial_angle)/Get_Periodicity())+1)/2.;
	}
	initial_angle=PI*(sin(2*PI*(Global::change_type[recurrent].counter)/Get_Periodicity())+1)/12.;
	Position_Standard_Change(initial_angle);
	
	Calculate_Global_Optima();
	
	Global::change_type[recurrent].counter++;
}
void Rotation_DBG::Chaotic_Change(){

	for(int i=0;i<num_peakorfun;i++){
		height[i]=Chaotic_Value(height[i],min_height,max_height);
		weight[i]=Chaotic_Value(weight[i],min_width,max_width);
	}
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[chaotic].counter++;
}
void Rotation_DBG::Small_Step_Change(){
	Height_Standard_Change();
	Width_Standard_Change();
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[small_step].counter++;
}
void Rotation_DBG::Large_Step_Change(){
	Height_Standard_Change();
	Width_Standard_Change();
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[large_step].counter++;
}
void Rotation_DBG::Recurrent_Noisy_Change(){

	double initial_angle;
	double height_range=max_height-min_height;
	double width_range=max_width-min_width;
	double noisy;
	for(int i=0;i<num_peakorfun;i++){
		initial_angle=(double)Get_Periodicity()*i/num_peakorfun;
		height[i]=Sin_Value_Noisy(Global::change_type[recurrent_noisy].counter,min_height,max_height,height_range,initial_angle,recurrent_noisy_servity);	
		weight[i]=Sin_Value_Noisy(Global::change_type[recurrent_noisy].counter,min_width,max_width,width_range,initial_angle,recurrent_noisy_servity);		
	}
	initial_angle=PI*(sin(2*PI*(Global::change_type[recurrent_noisy].counter)/Get_Periodicity())+1)/12.;
	noisy=recurrent_noisy_servity*Global::normal.Next();
	Position_Standard_Change(initial_angle+noisy);
	
	Calculate_Global_Optima();
	Global::change_type[recurrent_noisy].counter++;
}
double Rotation_DBG::Evaluation(const double *const x){
	Copy(genes,x,dimension);
	for(int i=0;i<num_peakorfun;i++){
		fit[i]=0;
		for(int j=0;j<dimension;j++)
			fit[i]+=(genes[j]-position[i][j])*(genes[j]-position[i][j]);
		fit[i]=sqrt(fit[i]/dimension);
		fit[i]=height[i]/(1+weight[i]*fit[i]);
	}
	return  Extremum(fit,num_peakorfun,MAX);
}

void Rotation_DBG::Width_Standard_Change(){
	double step;
	for(int i=0;i<num_peakorfun;i++){
		
		step=width_severity*Standard_Change(Get_Change_Type(),min_width,max_width);
		weight[i]=weight[i]+step;
		
		if(weight[i]>max_width||weight[i]<min_width) weight[i]=weight[i]-step;

	}
}
bool Rotation_DBG::Set_Periodicity(const int p){
	if(p<1) return false;
	Real_DBG::Set_Periodicity(p);
	return true;
}
void Rotation_DBG::Set_Width_Severity(const double sw){
	width_severity=sw;
}
void Rotation_DBG::Calculate_Global_Optima(){
	global_optima=Extremum(height,num_peakorfun,MAX);
	  for(int i=0;i<num_peakorfun;i++)
		 if(height[i]==global_optima) Copy(global_optima_position,position[i],dimension);
}
void Rotation_DBG::Dimension_Decrease(){
	Rotation_DBG* r_dbg=new Rotation_DBG(num_peakorfun,Global::num_dim);
	if(Global::change==recurrent||Global::change==recurrent_noisy)
		r_dbg->Set_Periodicity(periodicity);
	r_dbg->Parameter_Setting(*this);
	r_dbg->Calculate_Global_Optima();
	Global::g_dbg=r_dbg;					
}
void Rotation_DBG::Dimension_Increase(){
	Rotation_DBG* r_dbg=new Rotation_DBG(num_peakorfun,Global::num_dim);
	if(Global::change==recurrent||Global::change==recurrent_noisy)
		r_dbg->Set_Periodicity(periodicity);
	r_dbg->Parameter_Setting(*this);
	double lower,upper;
	lower=r_dbg->boundary[dimension].lower=Global::boundary.lower;
	upper=r_dbg->boundary[dimension].upper=Global::boundary.upper;

	for(int i=0;i<num_peakorfun;i++){
		r_dbg->position[i][dimension]= lower+(upper-lower)*Global::uniform.Next();						
		r_dbg->initial_position[i][dimension]=r_dbg->position[i][dimension];									
	}			
							
	if(Global::change==recurrent||Global::change==recurrent_noisy){	
		for(int i=0;i<r_dbg->periodicity;i++){
			if(Global::change_type[Global::change].counter<=i) break;
			for(int j=0;j<num_peakorfun;j++)
				r_dbg->rotation_planes[i][j][dimension]=dimension;
		}				
	}
	r_dbg->Calculate_Global_Optima();
	Global::g_dbg=r_dbg;
}

void Rotation_DBG::Parameter_Setting(General_DBG &g_dbg){
	Real_DBG::Parameter_Setting(g_dbg);
	width_severity=static_cast<Rotation_DBG*>(&g_dbg)->width_severity;
}

