#pragma once

class CPoint{
public:
	double *x;
	double *v;
public:
	CPoint();
	CPoint(const CPoint &p);
	CPoint &operator =(const CPoint & p);
	~CPoint();
	CPoint &operator =(const double * ve);
	double Distance(const CPoint &);
};

