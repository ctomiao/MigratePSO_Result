// Matrix.h
#pragma once
#include <fstream>
using namespace std;

class Matrix{								// *****************Orthogonal rotation matrix***********************
public:
	int cow,row;							// matrix size 
private:
	double **data;							// value of each element
public:
	Matrix();
	Matrix(const int c,const int r);
	~Matrix();
	Matrix & operator *(const Matrix & m);
	Matrix & operator=(const Matrix & m);
	bool Identity();
	void Set_Rotation(const int &r,const int &c,const double &angle);
	const double **Get_Data()const {
		return const_cast<const double **>(data);
	}
	void Set_Data(const double *d, const int &c,const int & r=1);
	// used for debug
	void Print(ofstream & out){
		for(int i=0;i<row;i++){
			for(int j=0;j<cow;j++)
				out<<data[i][j]<<" ";
			out<<endl;
	}
	};
};

