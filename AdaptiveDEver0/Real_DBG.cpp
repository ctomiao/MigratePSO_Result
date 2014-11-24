//Real_DBG.cpp
#include "stdafx.h"
#include "Real_DBG.h"

using namespace std;

Real_DBG::Real_DBG(const int num,const int dim):General_DBG(dim),num_peakorfun(num),min_height(Global::min_height),max_height(Global::max_height){
	genes= new double [dim];
	position= new double*[num_peakorfun];
	initial_position=new double*[num_peakorfun];
	for(int i=0;i<num_peakorfun;i++){
		position[i]=new double[dim];
		initial_position[i]=new double[dim];
	}
	height=new double[num_peakorfun];
	boundary=new Boundary[dim];
	height_severity=0;
	
	prediction=false;
	fit=new double [num_peakorfun];
	weight=new double[num_peakorfun];
	rotation_matrix= new Matrix [num_peakorfun];
	global_optima_position=new double[dim];
}
Real_DBG::~Real_DBG(){
	delete [] boundary;
	delete [] height;
	for(int i=0;i<num_peakorfun;i++){
		delete [] position[i];
		delete [] initial_position[i];
	}
	delete [] position;
	delete [] initial_position;

	delete [] genes;
	delete [] fit;
	delete [] weight;
	delete [] rotation_matrix;
	if(Get_Periodicity()>0){
		for(int j=0;j<Get_Periodicity();j++){
			for(int i=0;i<num_peakorfun;i++)
				delete [] rotation_planes[j][i];
			delete []rotation_planes[j];
		}
		delete []rotation_planes;
	}
	delete [] global_optima_position;
}
Real_DBG & Real_DBG::operator =(const Real_DBG &r_dbg){
	if(this==&r_dbg) return *this;
	
	General_DBG::operator =(r_dbg);
	num_peakorfun=r_dbg.num_peakorfun;
	Set_Boundary(r_dbg.boundary );
	Set_Height(r_dbg.height);

	Copy(genes,r_dbg.genes,dimension);
	for(int i=0;i<num_peakorfun;i++){
		Copy(position[i],r_dbg.position[i],dimension);
		Copy(initial_position[i],r_dbg.initial_position[i],dimension);
	}	
	height_severity=r_dbg.height_severity;
	prediction=r_dbg.prediction;
	Copy(fit,r_dbg.fit,num_peakorfun);
	Copy(weight,r_dbg.weight,num_peakorfun);
	
	Copy(rotation_matrix,r_dbg.rotation_matrix,num_peakorfun);				
	global_optima=r_dbg.global_optima;
	
	for(int i=0;i<periodicity;i++){
		for(int j=0;j<num_peakorfun;j++)
			Copy(rotation_planes[i][j],r_dbg.rotation_planes[i][j],dimension);
	}
	Copy(global_optima_position,r_dbg.global_optima_position,dimension);		
	return *this;
}
void Real_DBG::Set_Boundary(const Boundary *b){
		Copy(boundary,b,dimension);	
}
void Real_DBG::Set_Boundary(const double *b){
	for(int j=0;j<dimension;j++){
			boundary[j].upper=fabs(b[j]);
			boundary[j].lower=-fabs(b[j]);
	}
}
void Real_DBG::Set_Boundary(const Boundary &b){
	for(int j=0;j<dimension;j++){
		boundary[j].upper=b.upper;
		boundary[j].lower=b.lower;
	}
}
void Real_DBG::Set_Boundary(const double &b){
	for(int j=0;j<dimension;j++){
			boundary[j].upper=fabs(b);
			boundary[j].lower=-fabs(b);
	}
}
bool Real_DBG::Set_Periodicity(const int p){
	if(p<1) return false;
	General_DBG::Set_Periodicity(p);
	rotation_planes=new int**[p];
	for(int i=0;i<p;i++){
		rotation_planes[i]=new int*[num_peakorfun];
		for(int j=0;j<num_peakorfun;j++)
			rotation_planes[i][j]=new int[dimension];
	}
	return true;
}

void Real_DBG::Set_Height_Severity(const double hs){
	height_severity=hs;
}
void Real_DBG::Correction(){
	for(int j=0;j<dimension;j++){
		if(genes[j]>boundary[j].upper)
			genes[j]=boundary[j].upper;
		else if(genes[j]<boundary[j].lower)
			genes[j]=boundary[j].lower;
	}
}

void Real_DBG::Height_Standard_Change(){
	double step;
	for(int i=0;i<num_peakorfun;i++){
		step=height_severity*Standard_Change(Get_Change_Type(),min_height,max_height);
		height[i]=height[i]+step;
		if(height[i]>Global::max_height||height[i]<Global::min_height) height[i]=height[i]-step;
		
	}
}
void Real_DBG::Position_Standard_Change(double angle){

	// for each basic function of dimension n(even number) , R=R(l1,l2)*R(l3,l4)*....*R(li-1,li), 0<=li<=n
	
	if(Get_Change_Type()==chaotic){
		for(int i=0;i<num_peakorfun;i++)
			for(int j=0;j<dimension;j++)
				position[i][j]=Chaotic_Value(position[i][j],boundary[j].lower,boundary[j].upper);
		return;
	}
	int * d=new int[dimension];
	Matrix I;
	for(int i=0;i<num_peakorfun;i++){
		if((Get_Change_Type()==recurrent||Get_Change_Type()==recurrent_noisy)&&Global::change_type[Get_Change_Type()].counter>=Get_Periodicity()){
			Copy(d,rotation_planes[Global::change_type[Get_Change_Type()].counter%Get_Periodicity()][i],dimension);
		}
		else{
			Initialize_RandomArray(d,dimension);
			if(Get_Change_Type()==recurrent||Get_Change_Type()==recurrent_noisy)
			Copy(rotation_planes[Global::change_type[Get_Change_Type()].counter][i],d,dimension);
		}

		if((Get_Change_Type()==recurrent||Get_Change_Type()==recurrent_noisy)&&Global::change_type[Get_Change_Type()].counter%Get_Periodicity()==0)
			Copy(position[i],initial_position[i],dimension);

		I.Identity();
		for(int j=0;j+1<dimension;j+=2){
			if(Get_Change_Type()==small_step||Get_Change_Type()==large_step||Get_Change_Type()==u_random)
				angle=Standard_Change(Get_Change_Type(), -PI,PI);
			I.Set_Rotation(d[j],d[j+1],angle);
			if(j==0) rotation_matrix[i]=I;
			else
				rotation_matrix[i]=rotation_matrix[i]*I;
		}
		Matrix m(dimension,1);
		m.Set_Data(position[i],dimension);
		m=m*rotation_matrix[i];
		Copy(genes,m.Get_Data()[0],dimension);
		Correction();
		Copy(position[i],genes,dimension);
	}
	delete [] d;
}
double Real_DBG::Get_Global_Optima(){
	return global_optima;
}
void Real_DBG::Parameter_Setting(General_DBG &g_dbg){
	General_DBG::Parameter_Setting(g_dbg);
	Real_DBG *r_dbg=static_cast<Real_DBG *>(&g_dbg);
	int dim=Global::num_dim<r_dbg->dimension?Global::num_dim:r_dbg->dimension;

	height_severity=r_dbg->height_severity;
	prediction=r_dbg->prediction;						
	global_optima=r_dbg->global_optima;
	Copy(weight,r_dbg->weight,num_peakorfun);
	Copy(height,r_dbg->height,num_peakorfun);

	for(int i=0;i<num_peakorfun;i++){
		Copy(position[i],r_dbg->position[i],dim);						
		Copy(initial_position[i],r_dbg->initial_position[i],dim);										
	}			
	Copy(boundary,r_dbg->boundary,dim);						
	Copy(global_optima_position,r_dbg->global_optima_position,dim);							
	if(Global::change==recurrent||Global::change==recurrent_noisy){	
		for(int i=0;i<r_dbg->periodicity;i++){
			if(Global::change_type[Global::change].counter<=i) break;
			for(int j=0;j<num_peakorfun;j++){
				if(dim==Global::num_dim){// the number of dimensions decreases
					for(int m=0,k=0;k<dim;k++,m++)
						if(r_dbg->rotation_planes[i][j][m]==dim) {k--;continue;}
						else
							rotation_planes[i][j][k]=r_dbg->rotation_planes[i][j][m];
					
				}else 
					Copy(rotation_planes[i][j],r_dbg->rotation_planes[i][j],dim);
				
			}
			
		}				
	}
}
