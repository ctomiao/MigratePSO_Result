#include "stdafx.h"
# include "Matrix.h"
#include "Global.h"

Matrix::Matrix(){
	row=Global::num_dim;
	cow=row=Global::num_dim;
	if(row==0) return;

	data=new double*[row];
	for(int i=0;i<row;i++)
		data[i]=new double[cow];
}
Matrix::Matrix(const int c, const int r):cow(c),row(r){
	data=new double*[row];
	for(int i=0;i<row;i++)
		data[i]=new double[cow];
	
}
void Matrix::Set_Data(const double *d,const int &c,const int &r){
	if(r!=row) return;
	for(int i=0;i<row;i++)
		for(int j=0;j<cow;j++)
			data[i][j]=d[j];
}
Matrix::~Matrix(){
	for(int i=0;i<row;i++)
		delete [] data[i];
	delete [] data;
}
Matrix & Matrix::operator *(const Matrix &m){

	if(cow!=m.row) return *this;
	Matrix r(m.cow,row);
	for(int i=0;i<row;i++){
		for(int j=0;j<m.cow;j++){
			r.data[i][j]=0;
			for(int k=0;k<cow;k++)
				r.data[i][j]+=data[i][k]*m.data[k][j];
		}
	}
	*this=r;
	return *this;
}
Matrix&Matrix::operator =(const Matrix &m){
	if(row!=m.row||cow!=m.cow) return *this;
	if(this==&m) return *this;
	for(int i=0;i<row;i++)
		for(int j=0;j<cow;j++)
			data[i][j]=m.data[i][j];
	return *this;
}
bool Matrix::Identity(){
	if(row!=cow) return false;
	for(int i=0;i<row;i++)
		for(int j=0;j<cow;j++)
			if(j!=i) data[i][j]=0.;
			else data[i][j]=1.;
	return true;
}

void Matrix::Set_Rotation(const int &r, const int &c,const double &angle){
	Identity();
	data[r][r]=cos(angle);
	data[r][c]=-sin(angle);
	data[c][r]=sin(angle);
	data[c][c]=cos(angle);
}
