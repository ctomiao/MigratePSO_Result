#pragma once
#include "Point.h"
#include <fstream>
using namespace std;
class CParticle{
public:
	CPoint pself;
	CPoint pnearest; // the closest point of other swarm or local optimum if it is within their radius
	double fitness;
	int flag;		// the parameter in the CLPSO to indicate how long have not move ahead
	int *F;			// the F in CLPSO
	double distance;  /// the distance of the particle to the whole swarm


public:
	static double *vmax;
	//static double c1,c2; //accelerators 
	//static double w;     // inertia weight    /// need to change here,because the c1 and c2 are changalbe
public:
	CParticle();
	CParticle(const CParticle &p);
	~CParticle();
	void Initialization();
	void Fun_Obj();
	CParticle & operator =(const CParticle &p);
	double Velocity();
	int Comparison(const CParticle &p);
	double Distance(const CParticle &p);
	void PSO_Move(const CParticle & lbest , const CParticle &gbest,const double,const double,const double);

};
