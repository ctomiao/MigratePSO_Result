#include "stdafx.h"
#include "Point.h"
#include <math.h>
#include "Global.h"
#include "Composition_DBG.h"
#include "Rotation_DBG.h"

CPoint::CPoint(){
	const int d=Global::g_dbg->Get_Dimension();
	x=new double[d];
	v=new double[d];
}
CPoint::~CPoint(){
	delete []x;
	x=0;
	delete []v;
	v=0;
}
CPoint & CPoint::operator =(const CPoint &p){
	int i;
	const int d=Global::g_dbg->Get_Dimension();
	for(i=0;i<d;i++){
		x[i]=p.x[i];
		v[i]=p.v[i];
	}
	return *this;
}
CPoint &CPoint::operator =(const double * ve){
	int i;
	const int d=Global::g_dbg->Get_Dimension();
	for(i=0;i<d;i++){
		v[i]=ve[i];
	}
	return *this;
}
double CPoint::Distance(const CPoint & position){
	double d=0;
	for(int i=0;i<Global::g_dbg->Get_Dimension();i++){
		d+=(x[i]-position.x[i])*(x[i]-position.x[i]);
	}
	return sqrt(d);
}