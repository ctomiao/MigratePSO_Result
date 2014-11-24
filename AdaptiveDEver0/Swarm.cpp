#include "stdafx.h"
#include "Swarm.h"

#include <iostream>
#include "Global.h"

#include "Composition_DBG.h"
#include "Rotation_DBG.h"

using namespace std;

CParticle *CSwarm::global_best;

int CSwarm::pop_num=0;
CSwarm *CSwarm::sub_swarm=0;
int CSwarm::max_popsize=40;
const int CSwarm::max_popnum=5;
const int CSwarm::largest_popsize=60;///// no matter how to migrate, the popsize can not exceed this value
const int CSwarm::smallest_popsize=20;
CSwarm::CSwarm(void){
	// pop= new CParticle[popsize];
	// pop_best = new CParticle[popsize];
	// PC =new double[popsize];
}
CSwarm::CSwarm(const CSwarm &s){
	int i;
	popsize=s.popsize ;
	pop= new CParticle[popsize];
	pop_best = new CParticle[popsize];
	PC =new double[popsize];
	w=s.w;
	c1=s.c1;
	c2=s.c2;
	for( i=0;i<popsize;i++){
		pop[i]=s.pop[i];
		pop_best[i]=s.pop_best[i];
		PC[i]=s.PC[i];
	}
	stage=s.stage;
	Best=s.Best;
	type= s.type;
	avg_fitness=s.avg_fitness;
	std_error=s.std_error;

}
CSwarm::CSwarm(const int pop_size){
	popsize=pop_size;
	pop= new CParticle[pop_size];
	pop_best = new CParticle[pop_size];
	PC=new double[pop_size];
	w=0.9;
	c1=c2=2.0;
	stage=0;
}
CSwarm &CSwarm::operator= (const CSwarm &s){
	if(this==&s)return *this;

	int i;
	popsize=s.popsize ;
	for( i=0;i<popsize;i++){
		pop[i]=s.pop[i];
		pop_best[i]=s.pop_best[i];
		PC[i]=s.PC[i];////
	}
	w=s.w;
	c1=s.c1;
	c2=s.c2;
	Best=s.Best;	
	type=s.type;
	stage=s.stage;
	avg_fitness=s.avg_fitness;
	std_error=s.std_error;
	return *this;
}
CSwarm::~CSwarm(void){
	if(pop){
		delete [] pop;
		pop=0;
	}
	if(pop_best) {
		delete [] pop_best;
		pop_best=0;
	}
	if(PC)
	{
		delete [] PC;
		PC=0;
	}
}

void CSwarm::Initial(const int pop_size){
	int i;
	popsize=pop_size;
	w=0.9;
	c1=c2=2.0;
	pop= new CParticle[pop_size];
	for( i=0;i<pop_size;i++) pop[i].Initialization();

	pop_best = new CParticle[pop_size];
	for( i=0;i<pop_size;i++)
		pop_best[i]=pop[i];

	PC= new double[popsize];
	cal_PC();    //// get the PC
	for( i=0;i<pop_size;i++) CL_Slection(i);
	stage=0;
	Best=pop[Find_Best()];

}
void CSwarm::Reinital(){
	int i;
	for( i=0;i<popsize;i++) pop[i].Initialization();
	for( i=0;i<popsize;i++)
		pop_best[i]=pop[i];
	w=0.9;
	c1=c2=2.0;

	cal_PC();    //// get the PC
	for( i=0;i<popsize;i++) CL_Slection(i);
	stage=0;
	Best=pop[Find_Best()];

}
const int CSwarm::Find_Best(void)const
{
	int index_best=0;
	for(int i=1;i<popsize;i++){
		if(pop[i].Comparison(pop[index_best])==1){
			index_best=i;
		}
	}
	return index_best;
}
void CSwarm::cal_PC()
{
	for(int i=0;i<popsize;i++)
		PC[i] = 0.05+0.45*(exp(10*i/(popsize-1))-1)/(exp(10)-1);
}
void CSwarm::Move(double best,int &fit_eva,double &r_value)
{
	switch(type)
	{
	case 0:				//GPSO	
		GMove(best,fit_eva,r_value);
		//CLMove(best,fit_eva,r_value);
		break;
	case 1:				/// LPSO
		LMove(best,fit_eva,r_value);
		//CLMove(best,fit_eva,r_value);
		break;
	case 2:				///APSO
		AMove(best,fit_eva,r_value);
		//CLMove(best,fit_eva,r_value);
		break;
	case 3:				///BBPSO
		BBMove(best,fit_eva,r_value);
		//CLMove(best,fit_eva,r_value);
		break;
	case 4:				/// CLPSO
		CLMove(best,fit_eva,r_value);
		//CLMove(best,fit_eva,r_value);
		break;
	default:
		break;
	}

}
void CSwarm::GMove(double best,int &fit_eva,double &r_value)
{
	c1=-2.0*Global::generation/Global::Max_gen+2.5;
	c2=2.0*Global::generation/Global::Max_gen+0.5;
	w=0.9-Global::generation*0.5/Global::Max_gen;
	for(int i=0;i<popsize;i++){
		pop[i].PSO_Move(pop_best[i],Best,w,c1,c2);
		record_info(best,fit_eva,r_value,i);
		Update_Best_Pbest();

	}
}
void CSwarm::LMove(double best,int &fit_eva,double &r_value)
{
	int i;
	int left,right,lbest;
	c1=-2.0*Global::generation/Global::Max_gen+2.5;
	c2=2.0*Global::generation/Global::Max_gen+0.5;
	w=0.9-Global::generation*0.5/Global::Max_gen;
	for(i=0;i<popsize;i++)
	{
		if(i==0)
		{
			left=1;
			right=popsize-1;
		}
		else
		{
			left=(i+1)%popsize;
			right=(i-1)%popsize;
		}
		if(pop_best[i].Comparison(pop_best[left])==-1 && pop_best[right].Comparison(pop_best[left])==-1)
			lbest=left;
		else if(pop_best[i].Comparison(pop_best[right])==-1 && pop_best[left].Comparison(pop_best[right])==-1)
			lbest=right;
		else
			lbest=i;
		pop[i].PSO_Move(pop_best[i],pop_best[lbest],w,c1,c2);
		record_info(best,fit_eva,r_value,i);
		Update_Best_Pbest();
	}
}
void CSwarm::AMove(double best,int &fit_eva,double &r_value)
{
	evolution_state_estimation();
	for(int i=0;i<popsize;i++)
	{
		pop[i].PSO_Move(pop_best[i],Best,w,c1,c2);
		record_info(best,fit_eva,r_value,i);
		Update_Best_Pbest();
	}
}
void CSwarm::BBMove(double best,int &fit_eva,double &r_value)
{
	int D=Global::g_dbg->Get_Dimension();
	int i,j,a,b,c;
	int left,right,lbest;
	double u,alpha;
	double mincoordinate,maxcoordinate;
	for(i=0;i<popsize;i++)
	{
		if(i == 0)
		{
			left =1;
			right=popsize-1;
		}
		else
		{
			left = (i+1)%popsize;
			right = (i-1)%popsize;
		}
		if(pop_best[i].Comparison(pop_best[left])==-1 && pop_best[right].Comparison(pop_best[left])==-1)
			lbest=left;
		else if(pop_best[i].Comparison(pop_best[right])==-1 && pop_best[left].Comparison(pop_best[right])==-1)
			lbest=right;
		else
			lbest=i;
		if(pop[i].Comparison(Best)==0)
		{
			for(j=0;j<D;j++)
			{
				mincoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
				maxcoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;
				if(Global::uniform.Next()<0.5)
				{
					do{a=popsize*Global::uniform.Next();}while(a==i);
					do{b=popsize*Global::uniform.Next();}while(b==i||a==b);
					do{c=popsize*Global::uniform.Next();}while(c==i||c==a||c==b);
					pop[i].pself.x[j]=pop_best[a].pself.x[j]+0.5*(pop_best[b].pself.x[j]-pop_best[c].pself.x[j]);
					if(pop[i].pself.x[j]>maxcoordinate||pop[i].pself.x[j]<mincoordinate) pop[i].pself.x[j]=mincoordinate+(maxcoordinate-mincoordinate)*Global::uniform.Next();
				}
				else
					pop[i].pself.x[j]=pop_best[lbest].pself.x[j];
			}
		}
		else
		{
			for(j=0;j<D;j++)
			{
				mincoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
				maxcoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;
				u=(Best.pself.x[j]+pop_best[i].pself.x[j])/2;
				alpha=fabs(Best.pself.x[j]-pop_best[i].pself.x[j]);
				if(alpha==0)
					alpha=0.001;
				pop[i].pself.x[j]=rnorm2(u,alpha);
				if(pop[i].pself.x[j]>maxcoordinate||pop[i].pself.x[j]<mincoordinate) pop[i].pself.x[j]=mincoordinate+(maxcoordinate-mincoordinate)*Global::uniform.Next();
			}
		}
		pop[i].fitness=static_cast<Real_DBG*>(Global::g_dbg)->Evaluation(pop[i].pself.x);
		record_info(best,fit_eva,r_value,i);
		Update_Best_Pbest();
	}

}
void CSwarm::CLMove(double best,int &fit_eva,double &r_value)
{
	int i,j;
	int D=Global::g_dbg->Get_Dimension();
	double mincoordinate;
	double maxcoordinate;
	w=0.9-Global::generation*0.5/Global::Max_gen;
	for(i=0;i<popsize;i++)
	{
		if(pop[i].flag>=7)
		{

			CL_Slection(i);
			pop[i].flag=0;
		}
		for(j=0;j<D;j++)
		{
			mincoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
			maxcoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;
			pop[i].pself.v[j]=w*pop[i].pself.v[j]+
				1.49445*Global::uniform.Next()*(pop_best[pop[i].F[j]].pself.x[j]-pop[i].pself.x[j]);
			if(fabs(pop[i].pself.v[j])>CParticle::vmax[j]) 
				pop[i].pself.v[j]=-CParticle::vmax[j]+2*CParticle::vmax[j]*Global::uniform.Next();		
			pop[i].pself.x[j]=pop[i].pself.x[j]+pop[i].pself.v[j];
			if(pop[i].pself.x[j]>maxcoordinate||pop[i].pself.x[j]<mincoordinate)pop[i].pself.x[j]=mincoordinate+(maxcoordinate-mincoordinate)*Global::uniform.Next();
		}
		pop[i].fitness=static_cast<Real_DBG*>(Global::g_dbg)->Evaluation(pop[i].pself.x);
		record_info(best,fit_eva,r_value,i);
		Update_Best_Pbest();
	}
}
void CSwarm::jDEMove(double best,int &fit_eva,double &r_value)
{
}


void CSwarm::record_info(double best,int &fit_eva,double &r_value,int i)
{
	fit_eva++;
	if(fit_eva%(Global::sample_frequency*max_popnum)==0){
		if(Global::optimization_type==MIN)
			r_value+=(1-best/CSwarm::global_best->fitness);
		else
			r_value+=(1-CSwarm::global_best->fitness/best);
	}
	if(pop[i].Comparison(pop_best[i])==1)
	{
		pop_best[i]=pop[i];
		pop[i].flag=0;
	}
	else
		pop[i].flag++;
	// if((pop[i].flag>5|| pop[i].Distance(Best)==0.0001)&& pop[i].Comparison(Best)==-1)
	// {
	// pop[i].Initialization();
	// CL_Slection(i);
	// }

}
void CSwarm::Statistic(void)
{
	Calculate_Avgfit();
	std_error=0;
	for(int i=0;i<popsize;i++)
		std_error+=(pop[i].fitness-avg_fitness)*(pop[i].fitness-avg_fitness);
	std_error=sqrt(std_error/popsize);
}
void CSwarm::Calculate_Avgfit(void){
	avg_fitness=0;
	// calculate the average fitness of its swarm
	for(int i=0;i<popsize;i++)
		avg_fitness+=pop[i].fitness;
	avg_fitness/=popsize;
}
void CSwarm::evolution_state_estimation()
{
	int i,j,k;
	double d[largest_popsize][largest_popsize];    /// need to consider here
	double min_d,max_d;
	double F;
	int D=Global::g_dbg->Get_Dimension();
	for(i=0;i<popsize;i++){
		d[i][i] = 0;
		for(j=i+1;j<popsize;j++){
			d[i][j] = 0;
			for(k=0;k<D;k++){
				d[i][j]+=(pop[i].pself.x[k]-pop[j].pself.x[k])
					*(pop[i].pself.x[k]-pop[j].pself.x[k]);
			}
			d[i][j] = sqrt(d[i][j]);
			d[j][i] = d[i][j];
		}
	}
	min_d = max_d = 0;
	for(i=0;i<popsize;i++){
		pop[i].distance = 0;
		for(j=0;j<popsize;j++)
			pop[i].distance += d[i][j];
		pop[i].distance /= (popsize-1); //excluding itself
		if(min_d==0||min_d>pop[i].distance)
			min_d = pop[i].distance;
		if(max_d<pop[i].distance)
			max_d = pop[i].distance;
	}
	///////////////////////////////////
	if(max_d==min_d)
		F=1;
	else
		F = (pop[Find_Best()].distance-min_d)/(max_d-min_d);
	adaptive_parameters(F);
}
void CSwarm::adaptive_parameters(double F)
{
	w = 1/(1+1.5*exp(-2.6*F));
	if(F<=0.2){//convergence
		change_c1_c2(0.5,0.5);
		elitist_learning();
		cur_state=3;
	}else if(F<=0.3){//overlap
		if(last_state==3||last_state==4){//convergence
			change_c1_c2(0.5,0.5);
			elitist_learning();
			cur_state=3;
		}else if(last_state==2||last_state==1){//exploitation
			change_c1_c2(0.5,-0.5);
			cur_state=2;
		}
	}else if(F<=0.4){//exploitation
		change_c1_c2(0.5,-0.5);
		cur_state=2;
	}else if(F<=0.6){//overlap
		if(last_state==2||last_state==3){//exploitation
			change_c1_c2(0.5,-0.5);
			cur_state=2;
		}else if(last_state==1||last_state==4){//exploration
			change_c1_c2(1.0,-1.0);
			cur_state=1;
		}
	}else if(F<=0.7){//exploration
		change_c1_c2(1.0,-1.0);
		cur_state=1;
	}else if(F<=0.8){//overlap
		if(last_state==1||last_state==2){//exploration
			change_c1_c2(1.0,-1.0);
			cur_state=1;
		}else if(last_state==4||last_state==3){//jumping out
			change_c1_c2(-1.0,1.0);
			cur_state=4;
		}
	}else{//jumping out
		change_c1_c2(-1.0,1.0);
		cur_state=4;
	}
	last_state=cur_state;
}
void CSwarm::change_c1_c2(double h1,double h2)
{
	double UP = 2.5; //The up bound of C1, C2
	double DOWN = 1.5;//The low bound of C1, C2
	double MAX = 4.0; //the up bound of c1+c2
	double MIN = 3.0; //the low bound of c1+c2
	double delta1 = 0.05; // the low bound of the random value in equation (11)
	double delta2 = 0.1; // the high bound of the random value in equation (11)
	double sum;

	c1 = c1+h1*(delta1+(delta2-delta1)*Global::uniform.Next());
	c2 = c2+h2*(delta1+(delta2-delta1)*Global::uniform.Next());
	if(c1<DOWN) c1 = DOWN; else if(c1>UP) c1 = UP;
	if(c2<DOWN) c2 = DOWN; else if(c2>UP) c2 = UP;
	sum = c1+c2;
	if(sum<MIN) {  c1=c1*MIN/sum; c2=c2*MIN/sum; }
	else if(sum>MAX){
		c1=c1*MAX/sum; c2=c2*MAX/sum;
	}
}
void CSwarm::elitist_learning()
{
	int j;
	//double value;
	CParticle pos;
	pos=Best;
	int D=Global::g_dbg->Get_Dimension();
	double Sita1 = 1.0; //the start sita in equation (13)
	double Sita2 = 0.1; // the end sita in equation (13)
	j=Global::uniform.Next()*D;
	double mincoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
	double maxcoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;
	pos.pself.x[j]=pos.pself.x[j]+(maxcoordinate-mincoordinate)*rnorm2(0,(Sita1-(Sita1-Sita2)*Global::generation/(Global::Max_gen-1)));
	//pos.x[j] = pos.x[j]+(uBound-lBound)*rnorm2(0,(Sita1-(Sita1-Sita2)*gen/(MAXGENS-1)));
	if(pos.pself.x[j]>maxcoordinate||pos.pself.x[j]<mincoordinate) pos.pself.x[j]=mincoordinate+(maxcoordinate-mincoordinate)*Global::uniform.Next();
	pos.fitness = static_cast<Real_DBG*>(Global::g_dbg)->Evaluation(pos.pself.x);
	if(pos.Comparison(Best)==1)
		Best=pos;
}

void CSwarm::CL_Slection(int index)
{
	int D=Global::g_dbg->Get_Dimension();
	int i,p,q;
	bool flag1=false;
	for(i=0;i<D;i++)
	{
		if(rnorm2(0,1)<PC[index])
		{
			flag1=true;
			do{p=(Global::uniform.Next())*popsize;}while(p==index);
			do{q=(Global::uniform.Next())*popsize;}while(q==index || p==q);
			if(pop_best[p].Comparison(pop_best[q])==1)
				pop[index].F[i]=p;
			else
				pop[index].F[i]=q;
		}
		else
			pop[index].F[i]=index;
	}
	if(!flag1)
	{
		do{p=(Global::uniform.Next())*popsize;}while(p==index);
		q=(Global::uniform.Next())*D;
		pop[index].F[q]=p;
	}
}
void CSwarm::Update_Best_Pbest()
{

	int b=Find_Best();
	if(pop[b].Comparison(Best)==1)
	{
		Best=pop[b];
		//if(stage>0) stage--;
		stage=0;

	}
	else
		stage++;
	// update the best one of all swarms
	if(Best.Comparison(*global_best)==1)
		global_best=&Best;
}
void CSwarm::Add_Pop(CSwarm & w){
	CSwarm *sw=new CSwarm[pop_num+1];
	int i=0;
	for( i=0;i<pop_num;i++){
		sw[i].Initial(sub_swarm[i].popsize);
		sw[i]=sub_swarm[i];
	}
	sw[i].Initial(w.popsize);
	sw[i]=w;
	delete [] sub_swarm;
	pop_num++;
	sub_swarm= sw;

}
void CSwarm::Delete_Pop(int index){
	if(pop_num==1) {
		delete  [] sub_swarm;
		sub_swarm=0;
		pop_num--;
		return ;
	}
	CSwarm *sw=new CSwarm[pop_num-1];
	int i,j;
	for( i=0,j=0;j<pop_num;i++,j++){
		if(j==index){
			i--;
			continue;	
		}
		sw[i].Initial(sub_swarm[j].popsize);
		sw[i]=sub_swarm[j];
	}
	delete [] sub_swarm;
	pop_num--;
	sub_swarm=sw;
}
void CSwarm::Add_Particle(const CParticle &add){
	CParticle *sw=new CParticle[popsize+1];
	CParticle *swbest=new CParticle[popsize+1];
	double *PC2=new double [popsize+1];
	int i;
	for( i=0;i<popsize;i++){
		sw[i]=pop[i];
		swbest[i]=pop_best[i];
	}
	sw[i]=add;
	swbest[i]=add;

	delete [] pop;
	pop=0;
	delete [] pop_best;
	pop_best=0;
	delete [] PC;
	PC=0;

	pop=sw;
	pop_best=swbest;
	PC=PC2;

	popsize++;
	Best=pop[Find_Best()];
	cal_PC();
	for(i=0;i<popsize;i++)
	{
		CL_Slection(i);
	}

}
void CSwarm::Delete_Particle(int del){
	if(popsize==1) {
		delete  [] pop;
		delete []pop_best;
		delete []PC;
		pop_best=0;
		pop=0;
		PC=0;
		popsize--;
		return;
	}
	CParticle *sw=new CParticle[popsize-1];
	CParticle *swbest=new CParticle[popsize-1];
	double *PC2=new double[popsize-1];
	int i,j;
	for( i=0,j=0;j<popsize;i++,j++){
		if(j==del){
			i--;
			continue;	
		}
		sw[i]=pop[j];
		swbest[i]=pop_best[j];
	}
	delete [] pop;
	pop=0;
	delete []pop_best;
	pop_best=0;
	delete []PC;
	PC=0;

	pop=sw;
	pop_best=swbest;
	PC=PC2;

	popsize--;
	Best=pop[Find_Best()];
	cal_PC();
	for(i=0;i<popsize;i++)
	{
		CL_Slection(i);
	}
}
