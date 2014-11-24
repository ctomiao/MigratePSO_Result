// AdaptiveDE.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


/*==========================================================================
//  Implementation of Generalized Benchmark Generator for Dynamic Problems
//
//  See the details of GDBG in "Generalized Dynamic Benchmark Generator for CEC'2009 
//  Competition on Dynamic Optimization"
//  
//  The GDBG was tested by PSO, the source sodes was implemented by Changhe Li
//
//  If you have any questions about the codes, please contact 
//  Changhe Li at cl160@le.ac.uk
===========================================================================*/


#include "Global.h"
#include "Composition_DBG.h"
#include "Rotation_DBG.h"
#include "Swarm.h"
#include <direct.h>
#include <fstream>
#include <iostream>
using namespace std;

char resultdir[100]="five_PSO_40";
void Initial_Change_Counter()
{
	for(int i=0;i<Num_Change_Type;i++)
	{
		Global::change_type[i].type=(Change_type) i;
		Global::change_type[i].counter=0;
	}
}
void PSO_Initial()
{
	//CSwarm::max_popsize=50;

	//particle velocity initial
	if(CParticle::vmax!=0)
	{
		delete [] CParticle::vmax;
		CParticle::vmax=0;
	}

	CParticle::vmax=new double[Global::num_dim];
	for(int i=0;i<Global::num_dim;i++){
		CParticle::vmax[i]=(Global::boundary.upper-Global::boundary.lower)/2.;
	}
}
void System_Initial(int t,double seed,const int num_dim){

	static LGM_mixed urng(seed);
	Random::Set(urng);

	Global::min_dimension=5;
	Global::max_dimension=15;
	Global::num_dim=num_dim;
	//Global::Swarm_type[0]=GPSO;Global::Swarm_type[1]=LPSO;Global::Swarm_type[2]=APSO;
	//Global::Swarm_type[3]=BBPSO;Global::Swarm_type[4]=CLPSO;

	while(Global::num_dim<5||Global::num_dim>15){
		cout<<"please input dimensions within [5,15]:";
		cin>>Global::num_dim;
		getchar();
	}
	Global::boundary.Set_Boundary(-5,5);
	if(t==Num_Change_Type){
		Global::flag_dimension_change=true;
		Global::change=(Change_type)(u_random);
	}
	else{
		Global::flag_dimension_change=false;
		Global::change=(Change_type)(t);
	}
	if(t==recurrent||t==recurrent_noisy)
		Global::periodicity=12;
	else
		Global::periodicity=0;

	Initial_Change_Counter();

	PSO_Initial();

}
void Continuous_Setting(Real_DBG *p){
	p->Set_Boundary(Global::boundary);
	double *t=new double[Global::num_peakorfun];

	for(int i=0;i<Global::num_peakorfun;i++){
		if(p->Get_Change_Type()==chaotic)
			t[i]=Global::min_height+(Global::max_height-Global::min_height)*Global::uniform.Next();
		else
			t[i]=50;
	}
	p->Set_Height(t);
	p->Set_Height_Severity(5);

	double **position;
	position=new double*[Global::num_peakorfun];
	for(int i=0;i<Global::num_peakorfun;i++)
		position[i]=new double[Global::num_dim];
	for(int i=0;i<Global::num_peakorfun;i++)
		for(int j=0;j<Global::num_dim;j++){
			position[i][j]=p->Get_Boundary()[j].lower+(p->Get_Boundary()[j].upper-p->Get_Boundary()[j].lower)*Global::uniform.Next();	
		}
		p->Set_Position(const_cast<const double **>(position));
		delete [] t;
		for(int i=0;i<Global::num_peakorfun;i++)
			delete []position[i];
		delete [] position;
}
void Rotation_DBG_Setting(Rotation_DBG *p){

	Continuous_Setting(p);
	p->Set_Width_Severity(0.5);
	p->Set_Weight(5);		// between (1,10)
	p->Calculate_Global_Optima();
}
void Composition_DBG_Setting(Composition_DBG * p,const int f){

	Continuous_Setting(p);

	Fun_Name *basic_fun=new Fun_Name[Global::num_peakorfun];
	switch(f){
	case 1:
		for(int i=0;i<Global::num_peakorfun;i++) basic_fun[i]=Sphere;
		break;
	case 2:
		for(int i=0;i<Global::num_peakorfun;i++) basic_fun[i]=Rastrigin;
		break;
	case 3:
		for(int i=0;i<Global::num_peakorfun;i++) basic_fun[i]=Griewank;
		break;
	case 4:
		for(int i=0;i<Global::num_peakorfun;i++) basic_fun[i]=Ackley;
		break;
	case 5:
		basic_fun[0]=Sphere;		basic_fun[1]=Sphere;
		basic_fun[2]=Rastrigin;		basic_fun[3]=Rastrigin;
		basic_fun[4]=Weierstrass;	basic_fun[5]=Weierstrass;
		basic_fun[6]=Griewank;		basic_fun[7]=Griewank;
		basic_fun[8]=Ackley;		basic_fun[9]=Ackley;
		break;
	}
	p->Set_Basic_Function(basic_fun);
	double *t=new double[Global::num_peakorfun];
	for(int i=0;i<Global::num_peakorfun;i++)t[i]=1.;
	p->Set_Coverge_Sevrity(t);
	p->Set_Stretch_Severity();		
	p->Set_Rotation_Matrix();
	p->Calculate_Global_Optima();
	delete []basic_fun;
	delete [] t;

}


void System_Setting(General_DBG * g_dbg, const int f){

	g_dbg->Set_Change_Frequency(Global::change_frequency);
	g_dbg->Set_Change_Type(Global::change);

	g_dbg->Set_Periodicity(Global::periodicity);

	g_dbg->Set_Dimension_Change(Global::flag_dimension_change);


	if(g_dbg->Get_Change_Type()==chaotic) {
		Global::chaotic_constant=3.67;
		while(Global::chaotic_constant>4||Global::chaotic_constant<1){
			cout<<"invalid value of chaotic_constant,reset please"<<endl;
			cin>>Global::chaotic_constant;
			getchar();
		}
	}

	if(g_dbg->Get_Change_Type()==recurrent_noisy)
		g_dbg->Set_RR_Severity(0.8f);

	if(Composition_DBG * p=dynamic_cast<Composition_DBG*>(g_dbg)){
		Composition_DBG_Setting(p,f);
	}else if(Rotation_DBG * p=dynamic_cast<Rotation_DBG*>(g_dbg)){
		Rotation_DBG_Setting(p);
	}
}
void Update_Global_Best(){
	CSwarm::global_best=&CSwarm::sub_swarm[0].Best;
	for(int k=1;k<CSwarm::pop_num;k++)
		if(CSwarm::global_best->Comparison(CSwarm::sub_swarm[k].Best)==-1)
			CSwarm::global_best=&CSwarm::sub_swarm[k].Best;
}
void Generate_Swarm(){

	for(int j=0;j<CSwarm::max_popnum;j++){

		CSwarm *newp= new CSwarm(CSwarm::max_popsize);
		newp->type=j;
		newp->Initial(newp->popsize);
		CSwarm::Add_Pop(*newp);
	}
	Update_Global_Best();
}
void Delete_Swarm(){
	for(int j=0;j<CSwarm::pop_num;j++){
		CSwarm::Delete_Pop(j);
		j--;
	}
}
void recalculate_swarm_size()
{
	int i,j,k;
	double pm;
	int migrate_in,migrate_out;
	int a;
	if(Global::generation!=0 && Global::generation%10==0)
	{
		for( k=0;k<CSwarm::pop_num;k++)						// migrate here
			CSwarm::sub_swarm[k].Statistic();
		//printf("the std_err of each swarm is:%g\t%g\t%g\t%g\t%g\n",CSwarm::sub_swarm[0].std_error,CSwarm::sub_swarm[1].std_error,CSwarm::sub_swarm[2].std_error,CSwarm::sub_swarm[3].std_error,CSwarm::sub_swarm[4].std_error);
		//printf("the stage of each swarm is:%d\t%d\t%d\t%d\t%d\n",CSwarm::sub_swarm[0].stage,CSwarm::sub_swarm[1].stage,CSwarm::sub_swarm[2].stage,CSwarm::sub_swarm[3].stage,CSwarm::sub_swarm[4].stage);
		pm=0.1+0.9*(exp((10.0*Global::generation)/(Global::Max_gen))-1)/(exp(10.0)-1);
		for(i=0;i<CSwarm::pop_num;i++)
			for(j=i+1;j<CSwarm::pop_num;j++)
			{
				if(Global::optimization_type==MIN)
				{
					if(CSwarm::sub_swarm[i].avg_fitness < CSwarm::sub_swarm[j].avg_fitness)
					{migrate_in=i;migrate_out=j;}
					else {migrate_in=j;migrate_out=i;}
				}
				if(Global::optimization_type==MAX)
				{
					if(CSwarm::sub_swarm[i].avg_fitness < CSwarm::sub_swarm[j].avg_fitness)
					{migrate_in=j;migrate_out=i;}
					else {migrate_in=i;migrate_out=j;}
				}
				if(Global::uniform.Next()<pm && CSwarm::sub_swarm[migrate_in].popsize<CSwarm::largest_popsize
					&&CSwarm::sub_swarm[migrate_out].popsize>CSwarm::smallest_popsize)
				{
					a=int(CSwarm::sub_swarm[migrate_out].popsize*Global::uniform.Next());/// choose a particle to migrate out randomly
					CSwarm::sub_swarm[migrate_in].Add_Particle(CSwarm::sub_swarm[migrate_out].pop[a]);
					CSwarm::sub_swarm[migrate_out].Delete_Particle(a);
				}
			}
	}
}
void recalculate_swarm_size2()    //// to make sure each swarm has the same number of particles
{	
	int i,j,k,migrate,migrate_in,migrate_out;
	for(i=0;i<CSwarm::pop_num;i++)
		for(j=i+1;j<CSwarm::pop_num;j++)
		{
			if(CSwarm::sub_swarm[i].popsize>CSwarm::max_popsize && CSwarm::sub_swarm[j].popsize<CSwarm::max_popsize)
			{
				migrate_out=CSwarm::sub_swarm[i].popsize-CSwarm::max_popsize;
				migrate_in=CSwarm::max_popsize-CSwarm::sub_swarm[j].popsize;
				if(migrate_out>migrate_in)
					migrate=migrate_in;
				else
					migrate=migrate_out;
				for(k=0;k<migrate;k++)
				{
					CSwarm::sub_swarm[j].Add_Particle(CSwarm::sub_swarm[i].pop[k]);
					CSwarm::sub_swarm[i].Delete_Particle(k);
				}
			}
			else if(CSwarm::sub_swarm[i].popsize<CSwarm::max_popsize && CSwarm::sub_swarm[j].popsize>CSwarm::max_popsize)
			{
				migrate_out=CSwarm::sub_swarm[j].popsize-CSwarm::max_popsize;
				migrate_in=CSwarm::max_popsize-CSwarm::sub_swarm[i].popsize;
				if(migrate_out>migrate_in)
					migrate=migrate_in;
				else
					migrate=migrate_out;
				for(k=0;k<migrate;k++)
				{
					CSwarm::sub_swarm[i].Add_Particle(CSwarm::sub_swarm[j].pop[k]);
					CSwarm::sub_swarm[j].Delete_Particle(k);
				}
			}
		}
}

void PSO(const int f,const int t,double ***particle_num,double **best_rel_fit,double **best_abs_fit,double **fit,double **relative,General_DBG * g_dbg,int *number_dimension,const int num_run){

	Generate_Swarm();

	//int max_gen;
	int fit_evas;
	double r_value;
	for(int i=0;i<Global::num_change;i++){
		Global::Max_gen=(Global::change_frequency*Global::num_dim/CSwarm::max_popsize)/SWARM_NUM;
		number_dimension[i]=Global::num_dim;
		if(num_run==0){// allocate memory to fit and relative for the first run
			fit[i]=new double[Global::Max_gen];
			relative[i]=new double[Global::Max_gen];
			particle_num[i]=new double *[Global::Max_gen];
			for(int j=0;j<Global::Max_gen;j++){
				particle_num[i][j]=new double [SWARM_NUM];
				fit[i][j]=0;
				relative[i][j]=0;
				for(int k=0;k<SWARM_NUM;k++)
					particle_num[i][j][k]=0;
			}
		}
		//int g;
		fit_evas=0;
		r_value=0;
		for(Global::generation=0;Global::generation<Global::Max_gen;Global::generation++){
			if(Global::generation % 50 == 0 )
			cout<<"Fun"<<f<<" t"<<t<<" "<<i<<" "<<Global::Max_gen<<" "<<Global::num_dim<<" "<<Global::generation<<" : "<<CSwarm::global_best->fitness<<" "<<static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima()<<endl;
			fit[i][Global::generation]+=CSwarm::global_best->fitness;
			if(Global::optimization_type==MIN)
				relative[i][Global::generation]+=static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima()/CSwarm::global_best->fitness;
			else
				relative[i][Global::generation]+=CSwarm::global_best->fitness/static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima();

			for(int k=0;k<CSwarm::pop_num;k++){
				particle_num[i][Global::generation][k]+=CSwarm::sub_swarm[k].popsize;
				CSwarm::sub_swarm[k].Move(static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima(),fit_evas,r_value);
				//printf("There are %d particles in %d and the best_fitness is:%g,the stage is:%d \n",CSwarm::sub_swarm[k].popsize,k,CSwarm::sub_swarm[k].Best.fitness,CSwarm::sub_swarm[k].stage);
				//printf("the stage of %d swarm is:%d\n",k,CSwarm::sub_swarm[k].stage);
			}
			recalculate_swarm_size();			/// migrate here
			for(int k=0;k<CSwarm::pop_num;k++)
			{
				//if(CSwarm::sub_swarm[k].stage>=CSwarm::sub_swarm[k].popsize*10)
				// if(CSwarm::sub_swarm[k].popsize==CSwarm::smallest_popsize)
				//   CSwarm::sub_swarm[k].Reinital();
			}

		}

		if(Global::optimization_type==MIN)
			best_rel_fit[num_run][i]=static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima()/CSwarm::global_best->fitness;
		else
			best_rel_fit[num_run][i]=CSwarm::global_best->fitness/static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima();

		best_abs_fit[num_run][i]=fabs(CSwarm::global_best->fitness-static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima());
		best_rel_fit[num_run][i]=best_rel_fit[num_run][i]/(1+r_value/(Global::change_frequency*Global::num_dim/(Global::sample_frequency)));

		g_dbg->Change();

		if(g_dbg->Get_Dimension_Change_Flag()==true) {
			g_dbg->Dimension_Change();
			if(Composition_DBG * p=dynamic_cast<Composition_DBG*>(g_dbg)){
				p->Delete_Composition_DBG();

			}else if(Rotation_DBG * p=dynamic_cast<Rotation_DBG*>(g_dbg)){
				p->Delete_Rotation_DBG();
			}
			g_dbg=Global::g_dbg;
		}

		// if(g_dbg->Get_Dimension_Change_Flag()==false)
		// recalculate_swarm_size2();   
		Delete_Swarm();
		PSO_Initial();
		Generate_Swarm();

		// if(g_dbg->Get_Dimension_Change_Flag()==false){ 
		// for(int k=0;k<CSwarm::pop_num;k++)
		// {
		// CSwarm::sub_swarm[k].Reinital();
		// }
		// Update_Global_Best();
		//}
		// else{
		// Delete_Swarm();
		// PSO_Initial();
		// Generate_Swarm();
		// }

	}
	Delete_Swarm();

}
double marking(const int f, const int t){
	double mark;
	if(f==0){
		if(t==6) mark=0.2*0.1*0.5;
		else mark=0.2*0.15*0.5;
	}else{
		if(t==6) mark=0.16*0.1;
		else mark=0.16*0.15;
	}
	return mark;
}
void Output_result(char * file, double ***particle_num,double ** fit, double ** relative, double **best_rel_fit,double **best_abs_fit,int *number_dimension,int num_peak,ofstream & perfor,const int f, const int t,double &total_mark){

	char name[50];
	strcpy(name,file);
	strcat(name,"_fit.txt");
	ofstream ofit(name);

	strcpy(name,file);
	strcat(name,"_relative.txt");
	ofstream orel(name);

	strcpy(name,file);
	strcat(name,"_statistic.txt");
	ofstream osta(name);

	strcpy(name,file);
	strcat(name,"_particle_num.txt");
	ofstream pnum(name);

	for(int i=0;i<Global::num_change;i++){
		for(int j=0;j<(Global::change_frequency*number_dimension[i]/CSwarm::max_popsize)/SWARM_NUM;j++)
			//for(int j=0;j<Global::Max_gen;j++)
		{
			ofit<<fit[i][j]/Num_Run<<endl;
			orel<<relative[i][j]/Num_Run<<endl;
			for(int k=0;k<SWARM_NUM;k++)
				pnum<<particle_num[i][j][k]/Num_Run<<"\t";
			pnum<<endl;
		}
		ofit << endl << endl;
		orel << endl << endl;
	}
	double min,max,std=0,avg_rel=0, avg_abs=0;
	double avg_min=0,avg_max=0;
	for(int i=0;i<Num_Run;i++){
		min=max=best_abs_fit[i][0]; 
		for(int j=0;j<Global::num_change;j++){
			if(min>best_abs_fit[i][j])
				min=best_abs_fit[i][j];
			if(max<best_abs_fit[i][j])
				max=best_abs_fit[i][j];
			avg_rel+=best_rel_fit[i][j];
			avg_abs+=best_abs_fit[i][j];
		}
		avg_min+=min;
		avg_max+=max;
	}
	avg_min/=Num_Run;
	avg_max/=Num_Run;
	avg_abs=avg_abs/(Num_Run*Global::num_change);
	avg_rel=avg_rel/(Num_Run*Global::num_change);

	for(int i=0;i<Num_Run;i++){
		for(int j=0;j<Global::num_change;j++)
			std+=(best_abs_fit[i][j]-avg_abs)*(best_abs_fit[i][j]-avg_abs);	
	}
	std=sqrt(std/(Num_Run*Global::num_change-1));
	osta<<"avg_min: "<<avg_min<<endl;
	osta<<"avg_max: "<<avg_max<<endl;
	osta<<"avg_mean: "<<avg_abs<<endl;
	osta<<"std: "<<std<<endl;

	double node_mark=marking(f,t);
	perfor<<f+1<<"\t"<<t+1<<"\t"<<num_peak<<"\t"<<avg_abs<<"\t"<<avg_rel<<"\t"<<node_mark<<"\t"<<node_mark*avg_rel<<endl<<endl;
	total_mark+=node_mark*avg_rel;

	ofit.close();
	orel.close();
	osta.close();
	pnum.close();

	for(int i=0;i<Global::num_change;i++)
	{
		for(int j=0;j<(Global::change_frequency*number_dimension[i]/CSwarm::max_popsize)/SWARM_NUM;j++)
			delete [] particle_num[i][j];
		//delete [] particle_num [i];

	}
	for(int i=0;i<Global::num_change;i++)
	{

		delete [] particle_num [i];
		delete []fit[i];
		delete []relative[i];
	}
	delete [] particle_num;
	delete [] fit;
	delete [] relative;

}
void Generate_file_name(char * file, const int f, const int t, const int d, const int m){
	// the name returned by file, f: function, t: change type, d: dimension, m: the number of peaks or basci functions
	char t1[50],t2[50];
	char resultdir2[100];
	strcpy(resultdir2,resultdir);
	strcat(resultdir2,"/");
	sprintf(t1,"%d",m);
	strcpy(t2,"_peak");
	strcat(t2,t1);
	sprintf(t1,"%d",d);
	strcat(t1,t2);
	strcpy(t2,"_D");
	strcat(t2,t1);
	sprintf(t1,"%d",t+1);
	strcat(t1,t2);
	strcpy(t2,"_T");
	strcat(t2,t1);
	sprintf(t1,"%d",f+1);
	strcat(t1,t2);
	strcpy(t2,"F");
	strcat(t2,t1);
	strcpy(file,t2);
	strcat(resultdir2,file);
	strcpy(file,resultdir2);

}
void run(const int f, const int t, const int num_dim, const int num_peak, double *seed,ofstream &perfor,double & total_mark){
	//run(function, change_type, number_dimensions, number_peaks, seed, output file, performance mark)
	int *number_dimension=new int[Global::num_change]; // save the number of dimensions of each change
	double **relative=new double*[Global::num_change]; 
	double **fit=new double* [Global::num_change];
	double **best_rel_fit=new double *[Num_Run];
	double **best_abs_fit=new double *[Num_Run];
	double ***particle_num=new double **[Global::num_change];
	for(int i=0;i<Num_Run;i++){
		best_rel_fit[i]=new double[Global::num_change];	
		best_abs_fit[i]=new double[Global::num_change];
	}

	if(f==0){
		for(int i=0;i<Num_Run;i++){
			System_Initial(t,seed[i],num_dim);
			Global::optimization_type=MAX;
			Rotation_DBG *rot_dbg;
			Global::num_peakorfun=num_peak;
			Global::max_width=10.;
			Global::min_width=1;
			Global::min_height=10;
			Global::max_height=100;
			rot_dbg= Rotation_DBG::Get_Rotation_DBG();
			System_Setting(rot_dbg,f);
			Global::g_dbg=rot_dbg;
			PSO(f,t,particle_num,best_rel_fit,best_abs_fit,fit,relative,rot_dbg,number_dimension,i); // call PSO algorithm
			rot_dbg->Delete_Rotation_DBG();
		}
		char file[50];
		Generate_file_name(file,f,t,num_dim,Global::num_peakorfun); // generate output file name
		Output_result(file,particle_num, fit, relative, best_rel_fit,best_abs_fit,number_dimension,num_peak,perfor,f,t,total_mark);

	}else{
		for(int i=0;i<Num_Run;i++){
			System_Initial(t,seed[i],num_dim);
			Global::optimization_type=MIN;
			Composition_DBG * com_dbg;
			Global::num_peakorfun=10;
			Global::hn_s=2000.;
			Global::min_height=10;
			Global::max_height=100;

			com_dbg=Composition_DBG::Get_Composition_DBG();
			System_Setting(com_dbg,f);
			Global::g_dbg=com_dbg;

			PSO(f,t,particle_num,best_rel_fit,best_abs_fit,fit,relative,com_dbg,number_dimension,i); // call PSO algorithm
			com_dbg->Delete_Composition_DBG();
		}
		char file[50];
		Generate_file_name(file,f,t,num_dim,Global::num_peakorfun);// generate output file name
		Output_result(file,particle_num, fit, relative, best_rel_fit,best_abs_fit,number_dimension,num_peak,perfor,f,t,total_mark);
	}			
	for(int i=0;i<Num_Run;i++){
		delete [] best_rel_fit[i];
		delete [] best_abs_fit[i];
	}

	delete [] best_rel_fit;
	delete [] best_abs_fit;

	delete [] number_dimension;
}
int main(){
	srand(17);
	double seed[Num_Run];
	for(int i=0;i<Num_Run;i++)
		seed[i]=(double)rand()/RAND_MAX;
	//int num_dim;
	int num_peak;
	char resultdir2[50];
	strcpy(resultdir2,resultdir);
	strcat(resultdir2,"/performance.txt");
	_mkdir(resultdir);
	//ofstream perfor("PSO_dynamic_result_testrestart/performance");
	ofstream perfor(resultdir2);
	double total_mark=0;
	perfor<<"problem\tchange_type\tnum_peak\tavg_abs\tavg_rel\tcase_mark score\n"<<endl;
	for(int f=0;f<6;f++){// 6 test problems
		for(int t=0;t<Num_Change_Type+1;t++){// 7 change types
			if(f==0){
				for(int p=0;p<2;p++){// 2 test number of peak: 10 and 50
					switch(p){
					case 0:
						num_peak=10;
						break;
					case 1:
						num_peak=50;
						break;
					}
					run(f,t,10,num_peak,seed,perfor,total_mark);
				}
			}else	run(f,t,10,10,seed,perfor,total_mark);

		}
	}
	delete [] CParticle::vmax;
	perfor<<"total mark (100*sum(score)): "<<100*total_mark<<endl;
	perfor.close();
	return 1;
}