#include "stdafx.h"
#include "Particle.h"
#include "Global.h"
#include <iostream>

#include "Composition_DBG.h"
#include "Rotation_DBG.h"

using namespace std;

double *CParticle::vmax=0;
//double CParticle::c1=1.4960;
//double CParticle::c2=1.4960;
//double CParticle::w=0.7298f;

CParticle::CParticle(){
	int D=Global::g_dbg->Get_Dimension();
	F=new int[D];
	flag=0;
	//pself=new CPoint();
	//pnearest=new CPoint();
}
CParticle::~CParticle(){
	if(F)
	{
		delete [] F;
		F=0;
	}
	//delete pself;
	//delete pnearest;
}
void CParticle::Initialization(){
	int i;
	double u,l;
	int D=Global::g_dbg->Get_Dimension();
	if(F!=0){
		delete [] F;
		F=0;
	}
	F= new int[D];
	flag=0;
	for( i=0;i<D;i++){
		u=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[i].upper;
		l=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[i].lower;
		pself.x[i]=l+(u-l)*Global::uniform.Next();
		if(u-l>2*vmax[i]) pself.v[i]=-vmax[i]+2*vmax[i]*Global::uniform.Next();
		else pself.v[i]=l+(u-l)*Global::uniform.Next();
	}
	Fun_Obj();

}
void CParticle::Fun_Obj(){
	fitness=static_cast<Real_DBG*>(Global::g_dbg)->Evaluation(pself.x);
}
CParticle & CParticle::operator=(const CParticle &p){
	if(this==&p) return *this;
	int d=Global::g_dbg->Get_Dimension();
	pself=p.pself;
	fitness=p.fitness;
	pnearest=p.pnearest;
	flag=p.flag;
	for(int i=0;i<d;i++)
		F[i]=p.F[i];
	return *this;
}
CParticle::CParticle(const CParticle &p)
{
	pself=p.pself;
	fitness=p.fitness;
	pnearest=p.pnearest;
	int D=Global::g_dbg->Get_Dimension();
	flag=p.flag;
	F=new int[D];
	for(int i=0;i<D;i++)
		F[i]=p.F[i];
}

double CParticle::Velocity(){
	int i;
	double ve=0;

	for( i=0;i<Global::g_dbg->Get_Dimension();i++)
		ve+=pself.v[i]*pself.v[i];
	if(ve==0.0) return 0;
	return sqrt(ve);
}
void CParticle::PSO_Move(const CParticle & lbest, const CParticle &gbest,const double w,
						 const double c1,const double c2){
							 double mincoordinate,maxcoordinate;
							 for(int j=0;j<Global::g_dbg->Get_Dimension();j++){
								 mincoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
								 maxcoordinate=static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;

								 pself.v[j]=w*pself.v[j]+c1*Global::uniform.Next()*(lbest.pself.x[j]-pself.x[j])+c2*Global::uniform.Next()*(gbest.pself.x[j]-pself.x[j]);
								 if(fabs(pself.v[j])>vmax[j]) 
									 pself.v[j]=-vmax[j]+2*vmax[j]*Global::uniform.Next();

								 pself.x[j]=pself.x[j]+pself.v[j];
								 if(pself.x[j]>maxcoordinate||pself.x[j]<mincoordinate) pself.x[j]=mincoordinate+(maxcoordinate-mincoordinate)*Global::uniform.Next(); 
							 }
							 Fun_Obj();
}
int CParticle::Comparison(const CParticle &p){
	int flag1;
	switch(Global::optimization_type){
	case MIN:
		if(fitness<p.fitness) flag1=1; 
		else if(fitness==p.fitness) flag1=0;
		else flag1=-1;
		break;
	case MAX:
		if(fitness<p.fitness) flag1=-1;
		else if(fitness==p.fitness) flag1=0;
		else flag1=1;
		break;
	}
	return flag1;
}
double CParticle::Distance(const CParticle &p)
{
	return pself.Distance(p.pself);
}