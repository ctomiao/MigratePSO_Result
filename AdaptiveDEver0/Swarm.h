#pragma once

#include "Particle.h"
#include "Global.h"
#include <fstream>
using namespace std;

class CSwarm
{
public:
	int popsize;  // swarm size  
	int type;
	double avg_fitness; // agerage fitness of whole swarm
	double std_error;   //standard deviation of whole swarm

	CParticle Best;  // the best one of swarm
	CParticle *pop;	 // the population
	CParticle *pop_best; //population of previous best position

	const static int max_popnum;
	const static int largest_popsize;
	const static int smallest_popsize;
	double w;
	double c1,c2;

	double *PC;

	int cur_state,last_state;
	int stage;		/// to indicate how long have not move ahead
public:

	static CParticle *global_best; // global best one found by PSO

	// the total number of swarms
	static int pop_num; 
	static CSwarm *sub_swarm;
	static int max_popsize;

public:
	CSwarm(void);
	CSwarm(const CSwarm &s);
	CSwarm(const int pop_size);
	~CSwarm(void);
	void Initial(const int pop_size);
	void Reinital();
	const int Find_Best(void)const;

	void Move(double best, int&, double &);
	void GMove(double best, int&, double &);
	void LMove(double best, int&, double &);
	void AMove(double best, int&, double &);
	void BBMove(double best, int&, double &);
	void CLMove(double best, int&, double &);
	void jDEMove(double best, int&, double &);


public:
	void Statistic(void);
	// function for the CLPSO
	void CL_Slection(int index);    /// to calculate the F for the index particle;
	void cal_PC(); 					/// to calculate the PC;

	///function for the APSO
	void evolution_state_estimation();
	void adaptive_parameters(double);
	void change_c1_c2(double,double);
	void elitist_learning();

public:
	void Calculate_Avgfit(void);
	void Update_Best_Pbest();
	static void Add_Pop(CSwarm &w);
	static void Delete_Pop(int index);
	void Add_Particle(const CParticle &add);
	void Delete_Particle(int del);
	void record_info(double best, int&, double &,int);
	CSwarm &operator= (const CSwarm &s);
};