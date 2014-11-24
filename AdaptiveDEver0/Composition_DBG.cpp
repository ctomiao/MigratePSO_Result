//Composition_DBG.cpp
#include "stdafx.h"
#include "Composition_DBG.h"
Composition_DBG* Composition_DBG::cr_dbg=0;

Composition_DBG::Composition_DBG(const int num, const int dim):Real_DBG(num,dim){

	com_boundary=new Boundary*[num_basic_fun];
	for(int i=0;i<num_basic_fun;i++)
		com_boundary[i]=new Boundary[dimension];
	Set_ComBoundary();
	height_normalize_severity=Global::hn_s;
	
	converge_severity =new double[num_peakorfun];						
	stretch_severity= new double[num_peakorfun];
	component_function=new Fun_Name[num_peakorfun];
	
}
Composition_DBG::~Composition_DBG(){

	delete [] stretch_severity;
	delete [] converge_severity;

	for(int i=0;i<num_basic_fun;i++)
		delete [] com_boundary[i];
	delete [] com_boundary;
	delete [] component_function;
	
}
void Composition_DBG::Set_ComBoundary(){
	//Sphere=0,Rastrigin,Weierstrass,Griewank,Ackley
		for(int j=0;j<dimension;j++){
			com_boundary[Sphere][j].upper=100;
			com_boundary[Sphere][j].lower=-100;
			com_boundary[Rastrigin][j].upper=5;
			com_boundary[Rastrigin][j].lower=-5;
			com_boundary[Weierstrass][j].upper=0.5;
			com_boundary[Weierstrass][j].lower=-0.5;
			com_boundary[Griewank][j].upper=100;
			com_boundary[Griewank][j].lower=-100;
			com_boundary[Ackley][j].upper=32;
			com_boundary[Ackley][j].lower=-32;
		}
}
Composition_DBG & Composition_DBG::operator=(const Composition_DBG &com){
	if(this==&com) return *this;
	Real_DBG::operator =(com);
	Set_ComBoundary();
	Copy(converge_severity,com.converge_severity,num_peakorfun) ;						
	Copy(stretch_severity,com.stretch_severity,num_peakorfun);																				
	height_normalize_severity=com.height_normalize_severity;					
	Copy(component_function,com.component_function,num_peakorfun);	

	return *this;
}
void Composition_DBG::Set_Rotation_Matrix(){
	// for each basic function of dimension n(even number), R=R(l1,l2)*R(l3,l4)*....*R(ln-1,ln), 0<=li<=n
	Matrix I;
	
	int * d=new int[dimension];
	Initialize_RandomArray(d,dimension);
	for(int i=0;i<num_peakorfun;i++){
		for(int j=0;j+1<dimension;j+=2){
			double angle=2*PI*Global::uniform.Next();				// random angle for rotation plane of d[j]-d[j+1] from d[j]th axis to d[j+1]th axis
			I.Set_Rotation(d[j],d[j+1],angle);
			if(j==0) rotation_matrix[i]=I;
			else
				rotation_matrix[i]=rotation_matrix[i]*I;
			I.Identity();
		}
	}
	delete [] d;
}
void Composition_DBG::Set_Coverge_Sevrity( const double *cs){
	Copy(converge_severity,cs,num_peakorfun);
}
void Composition_DBG::Set_Stretch_Severity(){
	for(int i=0;i<num_peakorfun;i++)
		stretch_severity[i]=converge_severity[i]*(boundary[0].upper-boundary[0].lower)/(com_boundary[(int)component_function[i]][0].upper-com_boundary[(int)component_function[i]][0].lower);
}
void Composition_DBG::Set_Basic_Function(const Fun_Name *bf){
	Copy(component_function,bf,num_peakorfun);
}

double Composition_DBG::Evaluation(const double *const x){
	Copy(genes,x,dimension);

	
	for(int i=0;i<num_peakorfun;i++){ // calculate weight for each function
		weight[i]=0;
		for(int j=0;j<dimension;j++)
			weight[i]+=(genes[j]-position[i][j])*(genes[j]-position[i][j]);
		weight[i]=exp(-sqrt(weight[i]/(2*dimension*converge_severity[i]*converge_severity[i])));
	}
	
	for(int i=0;i<num_peakorfun;i++){ // calculate objective value for each function
		
		for(int j=0;j<dimension;j++)	// calculate the objective value of tranformation function i
			genes[j]=(genes[j]-position[i][j])/stretch_severity[i];//((1+fabs(position[i][j]/boundary[j].upper))*
		Matrix m(dimension,1);
		m.Set_Data(genes,dimension); 
		
		m=m*rotation_matrix[i];
		
		Copy(genes,m.Get_Data()[0],dimension);
		Correction(component_function[i]);
		fit[i]=Select_Fun(component_function[i]);
		
		for(int j=0;j<dimension;j++){ // calculate the estimate max value of funciton i
			genes[j]=boundary[j].upper;
			genes[j]/=stretch_severity[i];
		}
		m.Set_Data(genes,dimension);
		m=m*rotation_matrix[i];
		Copy(genes,m.Get_Data()[0],dimension);
		Correction(component_function[i]);
		double fmax=Select_Fun(component_function[i]);
		if(fmax!=0)
		fit[i]=height_normalize_severity*fit[i]/fabs(fmax);
			
		Copy(genes,x,dimension);
	}
	double sumw=0,wmax;
	wmax=Extremum(weight,num_peakorfun,MAX);
	for(int i=0;i<num_peakorfun;i++)
		if(weight[i]!=wmax)
			weight[i]=weight[i]*(1-pow(wmax,10));
	for(int i=0;i<num_peakorfun;i++)
		sumw+=weight[i];
	for(int i=0;i<num_peakorfun;i++)
		weight[i]/=sumw;
	double obj=0;
	for(int i=0;i<num_peakorfun;i++)
		obj+=weight[i]*(fit[i]+height[i]);
	return obj;
}
inline double Composition_DBG::Select_Fun(const Fun_Name &f){
	double value;
	switch(f){
	case Sphere:
		value=F_Sphere();
		break;
	case Rastrigin:
		value=F_Rastrigin();
		break;
	case Weierstrass:
		value=F_Weierstrass();
		break;
	case Griewank:
		value=F_Griewank();
		break;
	case Ackley:
		value=F_Ackley();
		break;
	default:
		break;
	}
	return value;
}
inline double Composition_DBG::F_Ackley(){
	double fitness=0;
	double s1=0,s2=0;
	for(int i=0;i<dimension;i++){
		s1+=genes[i]*genes[i];
		s2+=cos(2*PI*genes[i]);
	}
	fitness=-20*exp(-0.2*sqrt(s1/dimension))-exp(s2/dimension)+20+E;	
	return fitness;
}
inline double Composition_DBG::F_Griewank(){
	double s1=0,s2=1;
	for(int i=0;i<dimension;i++){
		s1+=genes[i]*genes[i]/4000.;
		s2*=cos(genes[i]/sqrt((double) (i+1)));
	}
	return s1-s2+1.;
}
inline double Composition_DBG::F_Rastrigin(){
	double fit=0;
	for(int i=0;i<dimension;i++)
		fit=fit+genes[i]*genes[i]-10.*cos(2*PI*genes[i])+10.;
	return fit;
}
inline double Composition_DBG::F_Sphere(){
	double fit=0;
	for(int i=0;i<dimension;i++)
		fit+=genes[i]*genes[i];
	return fit;
}
inline double Composition_DBG::F_Weierstrass(){
	double a=0.5,b=3;
	int kmax=20;
	double fit=0,s=0;
	for(int i=0;i<dimension;i++)
		for(int k=0;k<=kmax;k++)
			fit+=pow(a,k)*cos(2*PI*pow(b,k)*(genes[i]+0.5));
	for(int k=0;k<=kmax;k++)
			s+=pow(a,k)*cos(2*PI*pow(b,k)*0.5);
	s=s*dimension;
	return fit-s;
}
void Composition_DBG::Correction(const Fun_Name &f){
	for(int j=0;j<dimension;j++){
		if(genes[j]>com_boundary[f][j].upper)  genes[j]=com_boundary[f][j].upper;
		else if(genes[j]<com_boundary[f][j].lower)  genes[j]=com_boundary[f][j].lower;
	}
}

void Composition_DBG::Random_Change(){
	//change the global minimum value of each function
	Height_Standard_Change();	
	//change the position of global optimum of each function randomly
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[u_random].counter++;
}
void Composition_DBG::Recurrent_Noisy_Change(){
	double initial_angle;
	double height_range=max_height-min_height;
	
	double noisy;
	for(int i=0;i<num_peakorfun;i++){
		initial_angle=(double)Get_Periodicity()*i/num_peakorfun;
		height[i]=Sin_Value_Noisy(Global::change_type[recurrent_noisy].counter,min_height,max_height,height_range,initial_angle,recurrent_noisy_servity);	
	}
	initial_angle=PI*(sin(2*PI*(Global::change_type[recurrent_noisy].counter)/Get_Periodicity())+1)/12.;
	noisy=recurrent_noisy_servity*Global::normal.Next();
	Position_Standard_Change(initial_angle+noisy);
	
	Calculate_Global_Optima();
	Global::change_type[recurrent_noisy].counter++;
}
void Composition_DBG::Recurrent_Change(){

	double initial_angle;
	double height_range=max_height-min_height;

	for(int i=0;i<num_peakorfun;i++){
		initial_angle=(double)Get_Periodicity()*i/num_peakorfun;
		height[i]=min_height+height_range*(sin(2*PI*(Global::change_type[recurrent].counter+initial_angle)/Get_Periodicity())+1)/2.;
	}
	initial_angle=PI*(sin(2*PI*Global::change_type[recurrent].counter/Get_Periodicity())+1)/12.;
	Position_Standard_Change(initial_angle);
	
	Calculate_Global_Optima();
	Global::change_type[recurrent].counter++;
}
void Composition_DBG::Small_Step_Change(){

	Height_Standard_Change();
	Position_Standard_Change(0);
	
	Calculate_Global_Optima();
	Global::change_type[small_step].counter++;
}
void Composition_DBG::Large_Step_Change(){
	Height_Standard_Change();
	Position_Standard_Change(0);
	Calculate_Global_Optima();
	Global::change_type[large_step].counter++;
}
void Composition_DBG::Chaotic_Change(){

	for(int i=0;i<num_peakorfun;i++)
		height[i]=Chaotic_Value(height[i],min_height,max_height);	
	
	Position_Standard_Change(0);
	
	Calculate_Global_Optima();
	Global::change_type[chaotic].counter++;
}
void Composition_DBG::Calculate_Global_Optima(){
	global_optima=Extremum(height,num_peakorfun,MIN);
	  for(int i=0;i<num_peakorfun;i++)
		 if(height[i]==global_optima) Copy(global_optima_position,position[i],dimension);
}
void Composition_DBG::Dimension_Decrease(){
	Composition_DBG* r_dbg=new Composition_DBG(num_peakorfun,Global::num_dim);
	if(Global::change==recurrent||Global::change==recurrent_noisy)
		r_dbg->Set_Periodicity(periodicity);
	r_dbg->Parameter_Setting(*this);
	r_dbg->Calculate_Global_Optima();
	Global::g_dbg=r_dbg;
}
void Composition_DBG::Dimension_Increase(){
	Composition_DBG* r_dbg=new Composition_DBG(num_peakorfun,Global::num_dim);
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
void Composition_DBG::Parameter_Setting(General_DBG &g_dbg){
	Real_DBG::Parameter_Setting(g_dbg);
	Composition_DBG* r_dbg=dynamic_cast<Composition_DBG*>(&g_dbg);

	Copy(component_function,r_dbg->component_function,num_peakorfun);
	Set_Rotation_Matrix();
	Copy(converge_severity,r_dbg->converge_severity,num_peakorfun);
	Copy(stretch_severity,r_dbg->stretch_severity,num_peakorfun);

}

